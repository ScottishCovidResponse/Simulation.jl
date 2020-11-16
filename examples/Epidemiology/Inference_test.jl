## Load packages
using SimulationData
using Simulation
using Simulation.Units
using Unitful
using Unitful.DefaultSymbols
using DifferentialEquations
using DiffEqSensitivity
using Random
using Distributions
using Turing
using DataFrames
using StatsPlots

# include("../../src/Epidemiology/Inference.jl")

## Simple SIRGrowth model

### Set model parameters
#=
beta_env = transmission from environmental reservoir
beta_force = transmission from airborne force of infection
sigma = recovery
virus_growth = rate of virus produced per infected
virus_decay = decay rate of environmental reservoir
mean_dispersal_dist = average dispersal distance of virus per each infected
=#
param = (beta_env = 2.0/day,
         beta_force = 3.0/day,
         sigma = 0.02/day,
         virus_growth = 1e-3/day,
         virus_decay = 1e-3/day,
         mean_dispersal_dist = 5.0km)

### Set run parameters:
#=
times = amount of time to simulate over
interval = how often to record output
timestep = timestep of simulation
=#
runparams = (times = 2years,
             interval = 1day,
             timestep = 1day)

### Set grid size & area
grid_size = (4,4)
area = 100.0km^2


###  Run simple model
# This set up runs in ~0.2 seconds per repeat.

@time abuns = SIR_wrapper(grid_size, area, param, runparams) #1st run: 15sec
#=
Outputs an abundance matrix of compartment by grid cell over time. Compartments for the SIR model are: Susceptible, Infected, Recovered, Dead
=#

### Using ODEs to infer β_env, β_force, and σ

pop_size = 500_000.0 # value matching Inference.jl #35
Y = abuns

@model bayes_sirGrowth(y) = begin
    l = length(y) # number of timepoints
    β_env ~ Uniform(0.0,1.0) # transmission from environmental reservoir
    β_force ~ Uniform(0.0,1.0) # transmission from airborne force of infection
    σ ~ Uniform(0.0,1.0) #recovery

    virus_growth = 1e-3/day # rate of virus produced per infected
    virus_decay = 1e-3/day # decay rate of environmental reservoir
    mean_dispersal_dist = 5.0km # average dispersal distance of virus per each infected
    # I0 = .. Inference.jl #72
    I = pop_size*β_env + pop_size*β_force
    t0 = [pop_size-I, I, σ*I, 0.0] # Susceptible, Infected, Recovered, Dead at time 0
    # Set up
    Ncells = prod(grid_size)

    β_env = convert(typeof(1e-3/day), β_env/day)
    β_force = convert(typeof(1e-3/day), β_force/day)
    σ = convert(typeof(1e-3/day), σ/day)

    # Set initial population sizes for all pathogen categories
    virus = 0
    abun_v = DataFrame([
        (name="Environment", initial=virus),
        (name="Force", initial=0),
    ])
    numvirus = nrow(abun_v)

    # Set initial population sizes for all human categories
    susceptible = 500_000 * Ncells
    abun_h = DataFrame([
        (name="Susceptible", type=Susceptible, initial=susceptible),
        (name="Infected", type=Infectious, initial=0),
        (name="Recovered", type=Removed, initial=0),
        (name="Dead", type=Removed, initial=0),
    ])
    numclasses = nrow(abun_h)

    # Set non-pathogen mediated transitions
    transitions = DataFrame([
        (from="Infected", to="Recovered", prob=σ),
    ])

    # Set simulation parameters & create transition matrices
    birth = fill(0.0/day, numclasses)
    death = fill(0.0/day, numclasses)
    param = (birth = birth, death = death, virus_growth = virus_growth, virus_decay = virus_decay, beta_env = β_env, beta_force = β_force)

    # Set up simple gridded environment
    epienv = simplehabitatAE(298.0K, grid_size, area, NoControl())

    # Dispersal kernels for virus and disease classes
    dispersal_dists = fill(mean_dispersal_dist, Ncells)
    kernel = GaussianKernel.(dispersal_dists, 1e-10)
    movement = EpiMovement(kernel)

    # Traits for match to environment (turned off currently through param choice, i.e. virus matches environment perfectly)
    traits = GaussTrait(fill(298.0K, numvirus), fill(0.1K, numvirus))
    epilist = EpiList(traits, abun_v, abun_h, movement, transitions, param)

    # Create epi system with all information
    rel = Gauss{eltype(epienv.habitat)}()
    epi = EpiSystem(epilist, epienv, rel)

    # Seed infected category at a single location
    human(epi.abundances)[2, 1] = 100 * Ncells

    # params = [β_env,
    #           β_force,
    #           σ,
    #           virus_growth,
    #           virus_decay,
    #           mean_dispersal_dist]
    params = [epi]
    tspan = (0.0,float(l))
    prob = ODEProblem(ODE_wrapper,
            t0,
            tspan,
            epi)
    sol = solve(prob,
                Tsit5(),
                saveat = 1.0)
    sol_C = Array(sol)[4,:] # Cumulative cases
    sol_X = sol_C[2:end] - sol_C[1:(end-1)]

    # for i in 1:l
    #   y[i] ~ Poisson(sol_X[i])
    # end
  end;


ode_nuts = sample(bayes_sirGrowth(Y),NUTS(0.65),10000);


describe(ode_nuts)


plot(ode_nuts)


posterior = DataFrame(ode_nuts)
