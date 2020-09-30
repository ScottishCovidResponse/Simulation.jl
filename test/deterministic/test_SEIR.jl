using SparseArrays
using Unitful
using Unitful.DefaultSymbols
using Simulation
using Simulation.Units
using DataFrames
using Test
include("seirs_det.jl")
# parameters
# note: this is a Named Tuple, which is immutable, but very convenient for typing
N = 10_000
timestep = 1.0
params = (
	# general parameters
	I0 = inv(100.0), # initial proportion of infected individuals
	T = 30.0,   # time to run simulation
	c = N,    # population size
	α = 1/7,    # incubation rate (1/period)
	β = 1.0,    # transmission coefficient
	γ = 1/7,    # recovery rate
	σ = 0.0,    # immunity waning rate
	μ = 0.0,    # mortality rate
	ν = 0.0,    # disease induced mortality rate
	κ = 0.0,  # movement rate
)

soln_det = seirs_det(params)
plot_det(soln_det)

const stochasticmode = false
rngtype = stochasticmode ? Random.MersenneTwister : MedianGenerator
# Set up simple gridded environment
grid = (1, 1)
area = 1.0km^2
epienv = simplehabitatAE(298.0K, grid, area, NoControl())

# Set initial population sizes for all pathogen categories
abun_v = DataFrame([
    (name="Environment", initial=0),
    (name="Force", initial=0),
])
numvirus = nrow(abun_v)

# Set initial population sizes for all human categories
N = 10_000
exposed = 100
susceptible = N - exposed
infected = 0
abun_h = DataFrame([
    (name="Susceptible", type=Susceptible, initial=susceptible),
    (name="Exposed", type=OtherDiseaseState, initial=exposed),
    (name="Infected", type=Infectious, initial=infected),
    (name="Recovered", type=Removed, initial=0),
    (name="Dead", type=Removed, initial=0),
])
numclasses = nrow(abun_h)

# Set non-pathogen mediated transitions
mu = 1/7days
sigma = 1/7days
transitions = DataFrame([
  (from="Exposed", to="Infected", prob=mu),
  (from="Infected", to="Recovered", prob=sigma),
])

# Set simulation parameters
birth = fill(0.0/day, numclasses)
death = fill(0.0/day, numclasses)
beta_force = 1.0/day
beta_env = 0.0/day
virus_growth = 1.0/day
virus_decay = 1.0/day
param = (birth = birth, death = death, virus_growth = virus_growth, virus_decay = virus_decay, beta_env = beta_env, beta_force = beta_force)

# Dispersal kernels for virus and disease classes
dispersal_dists = fill(0.5km, prod(grid))
kernel = GaussianKernel.(dispersal_dists, 1e-10)
movement = EpiMovement(kernel)

# Traits for match to environment (turned off currently through param choice, i.e. virus matches environment perfectly)
traits = GaussTrait(fill(298.0K, numvirus), fill(0.1K, numvirus))
epilist = EpiList(traits, abun_v, abun_h, movement, transitions, param)

# Create epi system with all information
rel = NoRelContinuous{eltype(epienv.habitat)}()
epi = EpiSystem(epilist, epienv, rel, rngtype = rngtype)

# Run simulation
times = 31days; interval = 1day; timestep = 1day
abuns = zeros(Int64, numclasses, prod(grid), div(times, interval) + 1)
@time simulate_record!(abuns, epi, times, interval, timestep; save=false, save_path=".")

soln_det.soln.t .+= 1
plot_det(soln_det)
plot_epidynamics!(epi, abuns, linestyle = :dash, colours = [:blue, :orange, :green, :purple, :white])
