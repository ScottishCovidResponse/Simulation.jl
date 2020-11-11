using Simulation
using SimulationData
using Unitful
using Unitful.DefaultSymbols
using Simulation.Units
using Simulation.ClimatePref
using StatsBase
using Distributions
using AxisArrays
using HTTP
using Random
using DataFrames
using Plots
using CSV
using JLD

const stochasticmode = false
const seed = hash(time()) # seed used for Random.jl and therefore rngs used in Simulation.jl

Random.seed!(seed)

function run_model(api::DataPipelineAPI, times::Unitful.Time, interval::Unitful.Time, timestep::Unitful.Time, file::String; do_plot::Bool = false, do_download::Bool = true, save::Bool = false, savepath::String = pwd(), include_pollution = true)
    # Download and read in population sizes for Scotland
    scotpop = parse_scottish_population(api)

    # Read number of age categories
    age_categories = size(scotpop, 3)

    # Set initial population sizes for all pathogen categories
    abun_v = DataFrame([
        (name="Environment", initial=0),
        (name="Force", initial=fill(0, age_categories)),
    ])
    numvirus = sum(length.(abun_v.initial))

    # Set population to initially have no individuals
    abun_h = DataFrame([
        (name="Susceptible", type=Susceptible, initial=fill(0, age_categories)),
        (name="Exposed", type=Exposed, initial=fill(0, age_categories)),
        (name="Asymptomatic", type=Infectious, initial=fill(0, age_categories)),
        (name="Presymptomatic", type=Infectious, initial=fill(0, age_categories)),
        (name="Symptomatic", type=Infectious, initial=fill(0, age_categories)),
        (name="Hospitalised", type=OtherDiseaseState, initial=fill(0, age_categories)),
        (name="Recovered", type=Removed, initial=fill(0, age_categories)),
        (name="Dead", type=Removed, initial=fill(0, age_categories)),
    ])
    numclasses = nrow(abun_h)
    numstates = sum(length.(abun_h.initial))

    # Set up simple gridded environment
    area = (AxisArrays.axes(scotpop, 1)[end] + AxisArrays.axes(scotpop, 1)[2] -
        2 * AxisArrays.axes(scotpop, 1)[1]) *
        (AxisArrays.axes(scotpop, 2)[end] + AxisArrays.axes(scotpop, 2)[2] -
        2 * AxisArrays.axes(scotpop, 2)[1]) * 1.0

    # Sum up age categories and turn into simple matrix
    total_pop = dropdims(sum(Float64.(scotpop), dims=3), dims=3)
    total_pop = AxisArray(total_pop, AxisArrays.axes(scotpop)[1], AxisArrays.axes(scotpop)[2])
    total_pop.data[total_pop .≈ 0.0] .= NaN
    # Shrink to smallest bounding box. The NaNs are inactive.

    pollution = parse_pollution(api)
    pollution = pollution[5513m .. 470513m, 531500m .. 1221500m, :]
    pm2_5 = GriddedPollution(pollution[pollutant = "pm2-5"])
    pm2_5.matrix = shrink_to_active(pm2_5.matrix, .!isnan.(total_pop))
    pm2_5.matrix .-= mean(pm2_5.matrix)

    total_pop = shrink_to_active(total_pop);

    # Prob of developing symptoms
    p_s = fill(read_estimate(
               api,
               "human/infection/SARS-CoV-2/symptom-probability",
               "symptom-probability"
           ), age_categories)

    param_tab = read_table(api, "prob_hosp_and_cfr/data_for_scotland", "cfr_byage")
    # Prob of hospitalisation
    p_h = param_tab.p_h[1:end-1] # remove HCW
    pushfirst!(p_h, p_h[1]) # extend age categories
    append!(p_h, fill(p_h[end], 2)) # extend age categories
    # Case fatality ratio
    cfr_home = param_tab.cfr[1:end-1]
    pushfirst!(cfr_home, cfr_home[1])
    append!(cfr_home, fill(cfr_home[end], 2))
    cfr_hospital = param_tab.p_d[1:end-1]
    pushfirst!(cfr_hospital, cfr_hospital[1])
    append!(cfr_hospital, fill(cfr_hospital[end], 2))

    @assert length(p_s) == length(p_h) == length(cfr_home)

    # Time exposed
    T_lat = days(read_estimate(
        api,
        "human/infection/SARS-CoV-2/latent-period",
        "latent-period"
    )Unitful.hr)

    # Time asymptomatic
    T_asym = days(read_estimate(
        api,
        "human/infection/SARS-CoV-2/asymptomatic-period",
        "asymptomatic-period"
    )Unitful.hr)
    @show T_asym

    # Time pre-symptomatic
    T_presym = 1.5days
    # Time symptomatic
    T_sym = days(read_estimate(
        api,
        "human/infection/SARS-CoV-2/infectious-duration",
        "infectious-duration"
    )Unitful.hr) - T_presym
    # Time in hospital
    T_hosp = read_estimate(
        api,
        "fixed-parameters/T_hos",
        "T_hos"
    )days
    # Time to recovery if symptomatic
    T_rec = read_estimate(
        api,
        "fixed-parameters/T_rec",
        "T_rec"
    )days

    # Exposed -> asymptomatic
    mu_1 = (1 .- p_s) .* 1/T_lat
    # Exposed -> Pre-symptomatic
    mu_2 = p_s .* 1/T_lat
    # Pre-symptomatic -> symptomatic
    mu_3 = fill(1 / T_presym, age_categories)
    # Symptomatic -> hospital
    hospitalisation = p_h .* 1/T_sym
    # Asymptomatic -> recovered
    sigma_1 = (1 .- p_s) .* 1/T_asym
    # Symptomatic -> recovered
    sigma_2 = (1 .- p_h) .* (1 .- cfr_home) .* 1/T_rec
    # Hospital -> recovered
    sigma_hospital = (1 .- cfr_hospital) .* 1/T_hosp
    # Symptomatic -> death
    death_home = cfr_home .* 2/T_hosp
    # Hospital -> death
    death_hospital = cfr_hospital .* 1/T_hosp

    transitions = DataFrame([
        (from="Exposed", to="Asymptomatic", prob=mu_1),
        (from="Exposed", to="Presymptomatic", prob=mu_2),
        (from="Presymptomatic", to="Symptomatic", prob=mu_3),
        (from="Symptomatic", to="Hospitalised", prob=hospitalisation),
        (from="Asymptomatic", to="Recovered", prob=sigma_1),
        (from="Symptomatic", to="Recovered", prob=sigma_2),
        (from="Hospitalised", to="Recovered", prob=sigma_hospital),
        (from="Symptomatic", to="Dead", prob=death_home),
        (from="Hospitalised", to="Dead", prob=death_hospital)
    ])

    # Set simulation parameters
    birth_rates = fill(0.0/day, numclasses, age_categories)
    death_rates = fill(0.0/day, numclasses, age_categories)
    birth_rates[:, 2:4] .= uconvert(day^-1, 1/20years)
    death_rates[1:end-1, :] .= uconvert(day^-1, 1/100years)
    virus_growth_asymp = fill(0.1/day, age_categories)
    virus_growth_presymp = fill(0.5/day, age_categories)
    virus_growth_symp = fill(1.0/day, age_categories)
    virus_decay = 1.0/3days
    beta_force = fill(2.0/day, age_categories)
    beta_env = fill(1.0/day, age_categories)
    age_mixing = fill(1.0, age_categories, age_categories)
    param = (birth = birth_rates, death = death_rates, virus_growth = [virus_growth_asymp virus_growth_presymp virus_growth_symp], virus_decay = virus_decay, beta_force = beta_force, beta_env = beta_env, age_mixing = age_mixing)

    if include_pollution
        epienv = simplehabitatAE(298.0K, size(total_pop), area, Lockdown(20days), pollution = pm2_5)
        param = (; param..., pollution_infectivity = 1.1/(μg * m^-3))
    else
        epienv = simplehabitatAE(298.0K, size(total_pop), area, Lockdown(20days), pollution = NoPollution())
    end

    movement_balance = (home = fill(0.5, numclasses * age_categories), work = fill(0.5, numclasses * age_categories))

    # Dispersal kernels for virus and disease classes
    dispersal_dists = fill(1.0km, length(total_pop))
    thresholds = fill(1e-3, length(total_pop))
    kernel = GaussianKernel.(dispersal_dists, thresholds)
    home = AlwaysMovement(kernel)

    # Import commuter data (for now, fake table)
    active_cells = findall(.!isnan.(total_pop[1:end]))
    from = active_cells
    to = sample(active_cells, weights(total_pop[active_cells]), length(active_cells))
    count = round.(total_pop[to]/10)
    home_to_work = DataFrame(from=from, to=to, count=count)
    work = Commuting(home_to_work)
    movement = EpiMovement(home, work)

    # Traits for match to environment (turned off currently through param choice, i.e. virus matches environment perfectly)
    traits = GaussTrait(fill(298.0K, numvirus), fill(0.1K, numvirus))
    epilist = EpiList(traits, abun_v, abun_h, movement, transitions, param, age_categories, movement_balance)
    rel = Gauss{eltype(epienv.habitat)}()

    # multiple dispatch in action
    rngtype = stochasticmode ? Random.MersenneTwister : MedianGenerator

    ismissing(file) && error("No file supplied")
    locs = Int64.(CSV.read(file)[!, :location])
    initial_infecteds = 100
    seed_fun = seedfile!

    scotpop = shrink_to_active(scotpop, 3)
    cum_inf = Vector{Vector{Float64}}(undef, length(locs))
    for i in eachindex(locs)
        seeding = Simulation.EpiSeedInf(initial_infecteds, [locs[i]], seed_fun)

        # Create epi system with all information
        @time epi = EpiSystem(epilist, epienv, rel, permutedims(scotpop, (:age, :grid_x, :grid_y)), seeding, UInt16(1), rngtype = rngtype)

        # Populate susceptibles according to actual population spread
        cat_idx = reshape(1:(numclasses * age_categories), age_categories, numclasses)
        N_cells = size(epi.abundances.matrix, 2)

        # Turn off work moves for <20s and >70s
        epi.epilist.human.home_balance[cat_idx[1:2, :]] .= 1.0
        epi.epilist.human.home_balance[cat_idx[7:10, :]] .= 1.0
        epi.epilist.human.work_balance[cat_idx[1:2, :]] .= 0.0
        epi.epilist.human.work_balance[cat_idx[7:10, :]] .= 0.0

        # Run simulation
        abuns = zeros(UInt16, size(epi.abundances.matrix, 1), sum(epi.epienv.active), floor(Int, times/timestep) + 1)
        @time simulate_record!(abuns, epi, times, interval, timestep, save = save, save_path = joinpath(savepath, "locations/location_$i"))
        cum_inf[i] = sum(Float64, abuns[cat_idx[:, 2], :, :], dims = (1, 3))[1, :, 1]
    end

    # Write to pipeline
    #write_array(api, "simulation-outputs", "final-abundances", DataPipelineArray(abuns))

    if do_plot
        # View summed SIR dynamics for whole area
        category_map = (
            "Susceptible" => cat_idx[:, 1],
            "Exposed" => cat_idx[:, 2],
            "Asymptomatic" => cat_idx[:, 3],
            "Presymptomatic" => cat_idx[:, 4],
            "Symptomatic" => cat_idx[:, 5],
            "Hospital" => cat_idx[:, 6],
            "Recovered" => cat_idx[:, 7],
            "Deaths" => cat_idx[:, 8],
        )
        display(plot_epidynamics(epi, abuns, category_map = category_map))
        display(plot_epiheatmaps(epi, abuns, steps = [30]))
    end
    return cum_inf
end

config = "data_config.yaml"
# download_data_registry(config)

times = 1months; interval = 1day; timestep = 1day
file = "Top_100_locs.csv"

# Pollution run
abuns_pollution = StandardAPI(config, "test_uri", "test_git_sha") do api
    run_model(api, times, interval, timestep, file)
end;
JLD.save("Abuns_pollution_mean.jld", "abuns", abuns_pollution)

# Normal run
abuns_normal = StandardAPI(config, "test_uri", "test_git_sha") do api
    run_model(api, times, interval, timestep, file, include_pollution = false)
end;
JLD.save("Abuns_normal_mean.jld", "abuns", abuns_normal)

abuns_pollution = JLD.load("Abuns_pollution.jld", "abuns")
abuns_normal = JLD.load("Abuns_normal.jld", "abuns")

numclasses = 8
age_categories = 10
cat_idx = reshape(1:(numclasses * age_categories), age_categories, numclasses)
category_map = (
    "Susceptible" => cat_idx[:, 1],
    "Exposed" => cat_idx[:, 2],
    "Asymptomatic" => cat_idx[:, 3],
    "Presymptomatic" => cat_idx[:, 4],
    "Symptomatic" => cat_idx[:, 5],
    "Hospital" => cat_idx[:, 6],
    "Recovered" => cat_idx[:, 7],
    "Deaths" => cat_idx[:, 8],
)

# Compare epi dynamics
display(plot_epidynamics(epi, abuns_pollution, category_map = category_map, layout = (@layout [a{0.4w} b{0.6w}]), subplot = 1, title = "Pollution", size = (1200, 800), legend = false, margin = 10*Plots.mm))
display(plot_epidynamics!(epi, abuns_normal, category_map = category_map, subplot = 2, title = "No pollution", legend = :outerright, right_margin = 15 * Plots.mm, legendtitlefonthalign = :left))

# Plot proportion exposed against pollution
cum_inf_poll = sum(Float64, abuns_pollution[cat_idx[:, 2], :, :], dims = 1)[1, :, :]
for i in 1:size(cum_inf, 2)
    cum_inf[:, i] ./= total_pop[1:end]
end
prop_inf = mean(cum_inf, dims = 2)[:, 1]
poll = ustrip.(epi.epienv.pollution.matrix[1:end])
display(scatter(poll, prop_inf[1:end], xlab = "Pollution (\\mu g m^{-3})", ylab = "Proportion exposed", legend = false, zcolor = poll, mc = :default_r, msc = :white, ma = 0.8, size = (1000, 800), margin = 10*Plots.mm, ylim = (0, 1)))

cum_inf_poll = sum(Float64, abuns_pollution[cat_idx[:, 2], :, :], dims = (1, 3))[1, :, 1]
cum_inf = sum(Float64, abuns_normal[cat_idx[:, 2], :, :], dims = (1, 3))[1, :, 1]
prop_inf = cum_inf_poll ./ (cum_inf .+ cum_inf_poll)
poll = ustrip.(epi.epienv.pollution.matrix[1:end])
display(scatter(poll, prop_inf, xlab = "Pollution (\\mu g m^{-3})", ylab = "Proportion exposed", legend = false, zcolor = poll, mc = :default_r, msc = :white, ma = 0.8, size = (1000, 800), margin = 10*Plots.mm))


prop_inf1 = cum_inf[:, 10]
prop_inf2 = cum_inf[:, 30]
poll = ustrip.(epi.epienv.pollution.matrix[1:end])
display(scatter(poll, prop_inf1, xlab = "Pollution (\\mu g m^{-3})", ylab = "Proportion exposed", legend = false, zcolor = poll, mc = :default_r, msc = :white, ma = 0.8, size = (1000, 800), margin = 10*Plots.mm, ylim = (0, 1), layout = 2, subplot = 1, title = "Day 10"))
display(scatter!(poll, prop_inf2, xlab = "Pollution (\\mu g m^{-3})", ylab = "Proportion exposed", legend = false, zcolor = poll, mc = :default_r, msc = :white, ma = 0.8, size = (1000, 800), margin = 10*Plots.mm, ylim = (0, 1), subplot = 2, title = "Day 30"))

display(plot_epiheatmaps(epi, abuns_pollution, steps = [30]))


total_pop_mat = Matrix{Float64}(total_pop)
total_pop_mat[isnan.(total_pop)] .= 0
rank_pop = sortperm(Matrix{Int64}(total_pop_mat)[1:end], rev = true)
rank_pop = rank_pop[(total_pop_mat[rank_pop] .> 0)]
find_100 = rank_pop[abs.(total_pop_mat[rank_pop] .- 100) .< 20]
histogram(poll[find_100], xlab = "Pollution (\\mu g m^{-3})", ylab = "Count", legend = false)

using Diversity
pop = map(1:10) do i
    scotpop[:, :, i][find_100]
end
pop = Array(hcat(pop...)')
rho_bars = norm_sub_rho(Metacommunity(pop), 1.0)
max_rho = findmax(rho_bars[:diversity])[2]

ranked_rep = map(1:length(find_100)) do i
    norm_meta_rho(Metacommunity([pop[:, i] pop[:, max_rho]]), 1.0)[:diversity]
end

using CSV
using DataFrames
top_100 = find_100[sortperm(ranked_rep, rev = true)][1:100]
top_100_dat = DataFrame(location = top_100)
CSV.write("Top_100_locs.csv", top_100_dat)


top_100 = CSV.read("Top_100_locs.csv")
poll = ustrip.(epi.epienv.pollution.matrix[epi.epienv.active])
ref_array = fill(NaN, size(epi.epienv.active))
ref_array[epi.epienv.active] .= collect(1:sum(epi.epienv.active))
grid_ids = Int64.(ref_array[top_100[!, :location]])
cum_inf = hcat(abuns_pollution...)
cum_norm = hcat(abuns_normal...)
enough_cases = findall(cum_inf .> 5)
grids = [enough_cases[i][1] for i in eachindex(enough_cases)]
total = total_pop[epi.epienv.active]
display(scatter(poll[grids], cum_inf[enough_cases] ./total[grids], xlab = "Pollution (\\mu g m^{-3})", ylab = "Proportion exposed", legend = false, zcolor = poll[grids], mc = :default_r, msc = :white, ma = 0.8, size = (1000, 800), margin = 10*Plots.mm))
Plots.pdf("Top_100_of_total_mean.pdf")


total = cum_inf[enough_cases] ./ (cum_norm[enough_cases] .+ cum_inf[enough_cases])
display(scatter(poll[grids], total, xlab = "Pollution (\\mu g m^{-3})", ylab = "Proportion exposed", legend = false, zcolor = poll[grids], mc = :default_r, msc = :black, ma = 0.8, size = (1000, 800), margin = 10*Plots.mm))
Plots.pdf("Top_100_of_normal.pdf")
