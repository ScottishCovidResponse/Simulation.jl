"""
    seednone!(epi::EpiSystem, controls::NoControl, timestep::Unitful.Time)

Function that does no seeding of infected individuals.
"""
function seednone!(epi::EpiSystem, controls::NoControl, timestep::Unitful.Time)
    return controls
end

"""
    seedrandpop!(epi::EpiSystem, controls::Lockdown, timestep::Unitful.Time)

Function that seeds infection randomly, weighted by population size.
"""
function seedrandpop!(epi::EpiSystem, controls::Lockdown, timestep::Unitful.Time)
    rng = epi.abundances.rngs[Threads.threadid()]
    if (epi.seeding.initial_infected > 0) && (controls.current_date < controls.lockdown_date)
        inf = rand(rng, Poisson(epi.seeding.initial_infected * timestep /controls.lockdown_date))
        sus_ids = epi.epilist.human.susceptible
        exp_ids = epi.epilist.human.exposed
        w = weights(@view human(epi.abundances)[sus_ids, epi.ordered_active[1:epi.seeding.initial_infected]])
        pos = rand(rng, Multinomial(inf, w/sum(w)))
        human(epi.abundances)[sus_ids, epi.ordered_active[1:epi.seeding.initial_infected]] .-= reshape(pos, length(sus_ids), length(epi.ordered_active[1:epi.seeding.initial_infected]))
        human(epi.abundances)[exp_ids, epi.ordered_active[1:epi.seeding.initial_infected]] .+= reshape(pos, length(sus_ids), length(epi.ordered_active[1:epi.seeding.initial_infected]))
        human(epi.abundances)[human(epi.abundances) .< 0] .= 0
    elseif controls.current_date == controls.lockdown_date
        @info "Lockdown initiated - $(sum(human(epi.abundances)[epi.epilist.human.exposed, :])) individuals infected"
    end
    return controls
end

"""
    seedfile!(epi::EpiSystem, controls::Lockdown, timestep::Unitful.Time)

Function that seeds infection at locations described in a local file.
"""
function seedfile!(epi::EpiSystem, controls::Lockdown, timestep::Unitful.Time)
    rng = epi.abundances.rngs[Threads.threadid()]
    locs = epi.seeding.locs
    initial_inf = epi.seeding.initial_infected
    if (epi.seeding.initial_infected > 0) && (controls.current_date < controls.lockdown_date)
        inf = rand(rng, Poisson(initial_inf * timestep /controls.lockdown_date))
        sus_ids = epi.epilist.human.susceptible
        exp_ids = epi.epilist.human.exposed
        pos = rand(rng, Multinomial(inf, length(sus_ids)))
        for i in eachindex(locs)
            for j in eachindex(pos)
                if (pos[j] > 0) && (human(epi.abundances)[sus_ids[j], locs[i]] > 0)
                    human(epi.abundances)[sus_ids[j], locs[i]] -= pos[j]
                    human(epi.abundances)[exp_ids[j], locs[i]] += pos[j]
                end
            end
        end
        println(sum(human(epi.abundances)[exp_ids, :]))
        # human(epi.abundances)[sus_ids, locs] .-= reshape(pos, length(sus_ids), length(locs))
        # human(epi.abundances)[exp_ids, locs] .+= reshape(pos, length(exp_ids), length(locs))
        # human(epi.abundances)[human(epi.abundances) .< 0] .= 0
    elseif controls.current_date == controls.lockdown_date
        @info "Lockdown initiated - $(sum(human(epi.abundances)[epi.epilist.human.exposed, :])) individuals infected"
    end
    return controls
end
