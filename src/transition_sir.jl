mutable struct ViralLoad{U <: Unitful.Units} <: AbstractStateTransition
    location::Int64
    prob::TimeUnitType{U}
end

mutable struct Exposure{U <: Unitful.Units} <: AbstractStateTransition
    species::Int64
    location::Int64
    destination::Int64
    force_prob::TimeUnitType{U}
    virus_prob::TimeUnitType{U}
end

function getprob(rule::Exposure)
    return rule.force_prob, rule.virus_prob
end

mutable struct Infection{U <: Unitful.Units} <: AbstractStateTransition
    species::Int64
    location::Int64
    prob::TimeUnitType{U}
end

mutable struct Recovery{U <: Unitful.Units} <: AbstractStateTransition
    species::Int64
    location::Int64
    prob::TimeUnitType{U}
end

mutable struct SEIR{U <: Unitful.Units} <: AbstractStateTransition
    exposure::Exposure{U}
    infection::Infection{U}
    recovery::Recovery{U}
end

mutable struct ForceProduce{U <: Unitful.Units} <: AbstractPlaceTransition
    species::Int64
    location::Int64
    prob::TimeUnitType{U}
end

function getprob(rule::ForceProduce)
    return rule.prob
end

mutable struct ForceDisperse <: AbstractPlaceTransition
    species::Int64
    location::Int64
end

mutable struct Force{U <: Unitful.Units}
    forceprod::ForceProduce{U}
    forcedisp::ForceDisperse
end

function _run_rule!(epi::EpiSystem, rule::Exposure{U}, timestep::Unitful.Time) where U <: Unitful.Units
    rng = epi.abundances.rngs[Threads.threadid()]
    spp = getspecies(rule)
    loc = getlocation(rule)
    params = epi.epilist.params
    if epi.epienv.active[loc]
        N = sum_pop(epi.abundances.matrix, loc)
        env_inf = virus(epi.abundances)[1, loc] /
            (N^params.freq_vs_density_env)
        force_inf = virus(epi.abundances)[2, loc] /
            (N^params.freq_vs_density_force)
        expprob = (getprob(rule)[1] * force_inf + getprob(rule)[2] * env_inf) * timestep
        newexpprob = 1.0 - exp(-expprob)
        exposures = rand(rng, Binomial(epi.abundances.matrix[spp, loc], newexpprob))
        human(epi.abundances)[spp, loc] -= exposures
        human(epi.abundances)[spp + 1, loc] += exposures
    end
end

function _run_rule!(epi::EpiSystem, rule::Infection{U}, timestep::Unitful.Time) where U <: Unitful.Units
    rng = epi.abundances.rngs[Threads.threadid()]
    spp = getspecies(rule)
    loc = getlocation(rule)
    if epi.epienv.active[loc]
        infprob = getprob(rule) * timestep
        newinfprob = 1.0 - exp(-infprob)
        infs = rand(rng, Binomial(epi.abundances.matrix[spp, loc], newinfprob))
        human(epi.abundances)[spp, loc] -= infs
        human(epi.abundances)[spp + 1, loc] += infs
    end
end

function _run_rule!(epi::EpiSystem, rule::Recovery{U}, timestep::Unitful.Time) where U <: Unitful.Units
    rng = epi.abundances.rngs[Threads.threadid()]
    spp = getspecies(rule)
    loc = getlocation(rule)
    if epi.epienv.active[loc]
        recprob = getprob(rule) * timestep
        newrecprob = 1.0 - exp(-recprob)
        recs = rand(rng, Binomial(epi.abundances.matrix[spp, loc], newrecprob))
        human(epi.abundances)[spp, loc] -= recs
        human(epi.abundances)[spp + 1, loc] += recs
    end
end


function run_rule!(epi::EpiSystem, rule::SEIR{U}, timestep::Unitful.Time) where U <: Unitful.Units
    run_rule!(epi, rule.exposure, timestep)
    run_rule!(epi, rule.infection, timestep)
    run_rule!(epi, rule.recovery, timestep)
end

function _run_rule!(epi::EpiSystem, rule::Force{U}, timestep::Unitful.Time) where U <: Unitful.Units
    run_rule!(epi, rule.forceprod, timestep)
    run_rule!(epi, rule.forcedisp, timestep)
end

function _run_rule!(epi::EpiSystem, rule::ForceProduce{U}, timestep::Unitful.Time) where U <: Unitful.Units
    rng = epi.abundances.rngs[Threads.threadid()]
    spp = getspecies(rule)
    loc = getlocation(rule)
    # Calculate effective rates
    birthrate = getprob(rule) * timestep * human(epi.abundances)[spp, loc]
    births = rand(rng, Poisson(birthrate))
    # Spread force of infection over space
    if !iszero(births)
        virusmove!(epi, spp, loc, epi.cache.virusmigration, births)
    end
end

function _run_rule!(epi::EpiSystem, rule::ViralLoad{U}, timestep::Unitful.Time) where U <: Unitful.Units
    rng = epi.abundances.rngs[Threads.threadid()]
    params = epi.epilist.params
    loc = getlocation(rule)
    traitmatch = traitfun(epi, loc, 1)
    deathrate = getprob(rule) * timestep * traitmatch^-1
    # Convert death rate into 0 - 1 probability
    deathprob = 1.0 - exp(-deathrate)

    # Calculate how much virus degrades in the environment
    deaths = rand(rng, Binomial(virus(epi.abundances)[1, loc], deathprob))
    # Force of infection on average around half of timestep in environment
    survivalprob = exp(-deathrate/2.0)

    # So this much force of infection survives in the environment
    env_virus = rand(rng, Binomial(virus(epi.abundances)[2, loc], survivalprob * params.env_virus_scale))

    # Now update virus in environment and force of infection
    virus(epi.abundances)[1, loc] += env_virus - deaths
end

function _run_rule!(epi::EpiSystem, rule::ForceDisperse)
    rng = epi.abundances.rngs[Threads.threadid()]
    spp = getspecies(rule)
    loc = getlocation(rule)
    dist = Poisson(epi.cache.virusmigration[spp, loc])
    epi.cache.virusmigration[spp, loc] = rand(rng, dist)
end

function run_rule!(epi::EpiSystem, rule::R, timestep::Unitful.Time) where R <: AbstractTransition
    if typeof(rule) <: Exposure
        _run_rule!(epi, rule, timestep)
    elseif typeof(rule) <: Infection
        _run_rule!(epi, rule, timestep)
    elseif typeof(rule) <: Recovery
        _run_rule!(epi, rule, timestep)
    elseif typeof(rule) <: ForceDisperse
        _run_rule!(epi, rule)
    elseif typeof(rule) <: ViralLoad
        _run_rule!(epi, rule, timestep)
    elseif typeof(rule) <: ForceProduce
        _run_rule!(epi, rule, timestep)
    end
end

function create_transition_list(epilist::EpiList, epienv::GridEpiEnv)
    params = epilist.params

    state_list = [SEIR(Exposure(1, loc, params.transition_force[2, 1], params.transition_virus[2, 1]), Infection(2, loc, params.transition[3, 2]), Recovery(3, loc, params.transition[4, 3])) for loc in eachindex(epienv.habitat.matrix)]

    virus_list = [ViralLoad(loc, params.virus_decay) for loc in eachindex(epienv.habitat.matrix)]

    state_list = [virus_list; state_list]

    place_list = [ForceDisperse(spp, loc) for spp in eachindex(epilist.human.names) for loc in eachindex(epienv.habitat.matrix)]
    force_list = [ForceProduce(spp, loc, params.virus_growth[spp]) for spp in eachindex(epilist.human.names) for loc in eachindex(epienv.habitat.matrix)]
    place_list = [force_list; place_list]

    return TransitionList(state_list, place_list)
end

function new_update!(epi::EpiSystem, timestep::Unitful.Time)

    Threads.@threads for pl in epi.transitions.place
        run_rule!(epi, pl, timestep)
    end

    Threads.@threads for st in epi.transitions.state
        run_rule!(epi, st, timestep)
    end

    virus(epi.abundances)[2, :] .= sum(epi.cache.virusmigration, dims = 1)[1, :]

    # Invalidate all caches for next update
    invalidatecaches!(epi)

end


function new_simulate!(epi::E, duration::Unitful.Time, timestep::Unitful.Time) where E <: AbstractEpiSystem
  times = length(0s:timestep:duration)
  for i in 1:times
    new_update!(epi, timestep)
  end
end

function new_simulate_record!(storage::AbstractArray, epi::E,
  times::Unitful.Time, interval::Unitful.Time,timestep::Unitful.Time) where E <: AbstractEpiSystem
  ustrip(mod(interval,timestep)) == 0.0 || error("Interval must be a multiple of timestep")
  record_seq = 0s:interval:times
  time_seq = 0s:timestep:times
  storage[:, :, 1] = epi.abundances.matrix
  counting = 1
  for i in 2:length(time_seq)
    new_update!(epi, timestep);
    if time_seq[i] in record_seq
      counting = counting + 1
      storage[:, :, counting] = epi.abundances.matrix
    end
  end
  storage
end
