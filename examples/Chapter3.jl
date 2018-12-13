#addprocs(20)
#addprocs(20)

using Diversity
using Simulation
using RCall
using Distributions
using StatsBase
using ProgressMeter
using DataStructures
using Unitful.DefaultSymbols
using MyUnitful
using DataFrames
using JLD
# Set up initial parameters for ecosystem
numSpecies = 50

# Set up how much energy each species consumes
energy_vec = SimpleRequirement(repmat([2.0], numSpecies))

# Set probabilities
birth = 0.6/year
death = 0.6/year
l = 1.0
s = 0.5
boost = 1000.0
timestep = 1month

# Collect model parameters together (in this order!!)
param = EqualPop(birth, death, l, s, boost)

minT = -10.0°C
maxT = 10.0°C

grid = (20, 20)
area = 4.0km^2
totalK = 100000.0
individuals=10000

# Create ecosystem
kernel = GaussianKernel(10.0km, numSpecies, 10e-10)
movement = BirthOnlyMovement(kernel)

opts = rand(Normal(10.0, 10.0), numSpecies) * °C
vars = rand(Uniform(0, 25/9), numSpecies) * °C
traits = GaussTrait(opts, vars)
native = fill(true, numSpecies)
abun = rand(Multinomial(individuals, numSpecies))
sppl = SpeciesList(numSpecies, traits, abun, energy_vec,
movement, param, native)
abenv = simplehabitatAE(0.0°C,grid, totalK, area)
rel = Gauss{typeof(minT)}()
eco = Ecosystem(sppl, abenv, rel)

burnin = 200month; interval = 1month; reps = 1
lensim = length(0month:interval:times)
storage = generate_storage(eco, lensim, reps)
for j in 1:reps
if (j != 1) repopulate!(eco) end
  simulate!(eco, burnin, timestep)
  resetrate!(eco, 0.01°C/month)
  thisstore = view(storage, :, :, :, j)
  simulate_record!(thisstore, eco, times, interval, timestep)
end



using Diversity
using Simulation
using RCall
using Distributions
using StatsBase
using ProgressMeter
using DataStructures
using Unitful.DefaultSymbols
using MyUnitful
using DataFrames
using JLD
# Set up initial parameters for ecosystem
numSpecies = 50

# Set up how much energy each species consumes
energy_vec = SimpleRequirement(repmat([2.0], numSpecies))

# Set probabilities
birth = 0.6/year
death = 0.6/year
l = 1.0
s = 0.5
boost = 1000.0
timestep = 1month

# Collect model parameters together (in this order!!)
param = EqualPop(birth, death, l, s, boost)

minT = -10.0°C
maxT = 10.0°C

grid = (20, 20)
area = 4.0km^2
totalK = 100000.0
individuals=10000

# Create ecosystem
kernel = GaussianKernel(10.0km, numSpecies, 10e-10)
movement = BirthOnlyMovement(kernel)

function runsim(times::Int64)
  temps = readtable(Pkg.dir("Simulation", "examples", "Tempdat.csv"))
  temps = temps[sample(1:nrow(temps), 50), :]
  prefs = [temps[:OrigValueStr] temps[:OrigValueStr]+1]
  conv = collect(linspace(-5, 20, 10))
  prefs = mapslices(x -> conv[x], prefs, 1)
  opts = vcat(mapslices(x -> rand(Uniform(x[1], x[2])), prefs, 2)...) * °C
  #vars = collect(linspace(1, 4, numSpecies))
  vars = rand(Uniform(0, 25/9), numSpecies) * °C
  traits = GaussTrait(opts, vars)
  native = fill(true, numSpecies)
  abun = rand(Multinomial(individuals, numSpecies))
  sppl = SpeciesList(numSpecies, traits, abun, energy_vec,
   movement, param, native)
  abenv = tempgradAE(minT, maxT, grid, totalK, area, 0.0°C/month)
  rel = Gauss{typeof(minT)}()
  eco = Ecosystem(sppl, abenv, rel)

  steady_state = 100month; burnin = 200month; interval = 1month; reps = 1
  lensim = length(0month:interval:times)
  lensteady = length(0month:interval:steady_state)

  # Run simulations 10 times
  storage_steady = generate_storage(eco, lensteady, reps)
  storage_change = generate_storage(eco, lensim, reps)
  storage_newsteady = generate_storage(eco, lensteady, reps)
  for j in 1:reps
    if (j != 1) repopulate!(eco) end
      simulate!(eco, burnin, timestep)
      thisstoresteady = view(storage_steady, :, :, :, j)
      simulate_record!(thisstoresteady, eco, steady_state, interval, timestep)
      resetrate!(eco, 0.01°C/month)
      thisstorechange = view(storage_change, :, :, :, j)
      simulate_record!(thisstorechange, eco, times, interval, timestep)
      resetrate!(eco, 0.0°C/month)
      thisstorenewsteady = view(storage_newsteady, :, :, :, j)
      simulate_record!(thisstorenewsteady, eco, steady_state, interval, timestep)
  end
  #save("macrorun_2steady.jld", "storage_change", storage_change,
    #"storage_steady", storage_steady, "storage_newsteady", storage_newsteady)
end

runsim(1000month)