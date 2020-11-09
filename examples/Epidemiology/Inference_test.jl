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
    # params = [β_env,
    #           β_force,
    #           σ,
    #           virus_growth,
    #           virus_decay,
    #           mean_dispersal_dist]
    params = (beta_env = β_env,
             beta_force = β_force,
             sigma = σ,
             virus_growth = virus_growth,
             virus_decay = virus_decay,
             mean_dispersal_dist = mean_dispersal_dist,
             grid_size = grid_size,
             area = area)
    tspan = (0.0,float(l))
    prob = ODEProblem(ODE_wrapper,
            t0,
            tspan,
            params)
    sol = solve(prob,
                Tsit5(),
                saveat = 1.0)
    # sol_C = Array(sol)[4,:] # Cumulative cases
    sol_X = sol_C[2:end] - sol_C[1:(end-1)]

    # for i in 1:l
    #   y[i] ~ Poisson(sol_X[i])
    # end
  end;


ode_nuts = sample(bayes_sirGrowth(Y),NUTS(0.65),10000);


describe(ode_nuts)


plot(ode_nuts)


posterior = DataFrame(ode_nuts)
