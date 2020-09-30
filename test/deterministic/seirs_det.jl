using DifferentialEquations
using LabelledArrays
using Plots

function seirs!(du, u, params, t)
	# extract variables
	S, E, I, R = u
	N = S + E + I + R

	# extract parameters
	c = params.c
	α = params.α
	β = params.β
	γ = params.γ
	σ = params.σ
	μ = params.μ
	ν = params.ν

	# rates
	birthS     = μ * c
	deathS     = μ * S
	deathE     = μ * E
	deathI     = (μ + ν) * I
	deathR     = μ * R
	infection  = β * S * I / N
	incubation = α * E
	recovery   = γ * I
	waning     = σ * R

	# Susceptibles
	du[1] = + birthS - deathS - infection + waning
	# Exposed
	du[2] = - deathE + infection - incubation
	# Infectives
	du[3] = - deathI + incubation - recovery
	# Recovered
	du[4] = - deathR + recovery - waning
end #fn


function seirs_det(params)
	c = params.c
	T = params.T
	I0 = params.I0
	# initial conditions
	E0 = I0 * c
	u0 = [c - E0, E0, 0.0, 0.0]

	tspan = (0.0, T)
	p = params

	# set up the ODE and solve it
	prob = ODEProblem(seirs!, u0, tspan, p)
	soln = solve(prob, alg_hints=[:stiff])

	println("Final deterministic values:")
	println("[S,E,I,R] = ", round.(soln[end], digits=1))

	Imax, tidx = findmax(soln(soln.t, idxs=3))
	Imax, tmax = round.((Imax, soln.t[tidx]), digits=1)
	println("Peak infection: $Imax at time $tmax")

	(params = p, soln = soln)
end


function plot_det(results)
	X = results.soln
	T = results.params.T

	plot(X, lw=3,
	     xlim=(0, T),
	     xlab="Time",
		 ylab="Population",
		 title="Deterministic SEIRS model",
		 label=["S" "E" "I" "R"])
end #fn
