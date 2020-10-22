using Statistics
using StatsBase
using Distributions
using ProgressMeter
using SparseArrays

include("Person.jl")
include("ibm_fns.jl")

function seirs_ibm(params)
	states = (:S, :E, :I, :R)
	# main parameters
	T = params.T
	num_recs = params.num_recs
	num_itns = params.num_itns
	num_groups = params.num_groups
	method = params.method
	N = round(Int, params.c) # this must be an integer

	# population data
	X = [Persons.Person() for _ in 1:N]

	# initial conditions
	I0 = params.I0
	sampler_person = sampler(DiscreteUniform(1,N))
	primary_cases = sort([rand(sampler_person) for _ in 1:(I0 * N)])

	# data recording
	num_vars = 4 # 1 each for [S,E,I,R]
	rec_width = T / (num_recs - 1)
	records = zeros(num_vars, num_recs, num_itns)

	# dt is constant for now, it might vary in the future
	dt = min(rec_width, params.dt)

	# Contact Matrix times ("Who Contacted Who")
	size(params.WCW) == (N, N) || error("WCW matrix miss-sized")
	params.WCW .= 0

	# set up a nice progress bar
	values = ceil(Int, num_itns * T / dt)
	progbar = Progress(values, dt=1, barglyphs=BarGlyphs("[=> ]"))

	# Iterations start here
	for itn in 1:num_itns
		# reset to starting conditions
		reset_popn!(X)   # population
		t = 0.0          # time
		params.WCW .= 0       # contact matrix
		t_next_rec = 0.0 # time to next record
		rec = 0          # record number

		# add primary cases
		index_case = 0
		for i in primary_cases
			expose!(X[i], index_case, t, params)
		end #for

		# Time Loop begins here
		while t ≤ T
			next!(progbar)

			while t ≥ t_next_rec && rec < num_recs
				rec += 1
				t_next_rec += rec_width
				for (i, state) in enumerate(states)
					records[i, rec, itn] = get_sum(X, state)
				end
			end #while

			# simulate contacts and note exposures
			make_contacts!(params.WCW, X, t, dt, params)

			# update time
			t += dt

			# check for E → I → R → S and save S and I
			update_statuses!(X, t)
			if params.dropzero
				dropzeros!(params.WCW)
				nz = nonzeros(params.WCW)
				nz .-= one(eltype(nz))
			else
				nz = nonzeros(params.WCW)
				nz .= max.(nz, one(eltype(nz)))
				nz .-= one(eltype(nz))
			end
		end # time loop

		# finish recording in case a leap is > rec_width
		# alternatively could cap dt at rec_width, I just prefer it this way
		while rec < num_recs
			for (i, state) in enumerate(states)
				records[i, rec, itn] = get_sum(X, state)
			end
		end
	end # iterations

	(params = params, soln = records)
end #fn


function plot_ibm(results)
	# extract parameters to get time
	T = results.params.T
	num_recs = results.params.num_recs
	t = range(1, T + 1, length=num_recs)

	# simple mean of variables
	soln = results.soln
	X = dropdims(mean(soln, dims=3), dims=3)'

	plot(t, X, lw=3,
	xlim=(0, T),
	xlab="Time",
	ylab="Population",
	title="Individual-based SEIRS model",
	label=["S" "E" "I" "R"])
end #fn
