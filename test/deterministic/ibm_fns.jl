function reset_popn!(X)
	for x in X
		x.status = :S
		x.t_enter_E = Inf
		x.t_enter_I = Inf
		x.t_enter_R = Inf
		x.t_enter_S = Inf
		x.t_test = Inf
		x.t_isolated = Inf
		x.infected_by = -1
		x.alerted = Inf
		x.alerted_by = -1
	end #for

	nothing
end #fn



"""
    expose!(S, I, tE, params)
"""
function expose!(S, I, tE, params,
				 samplers = (α = sampler(Exponential(inv(params.α))),
	γ = sampler(Exponential(inv(params.γ))),
	σ = sampler(Exponential(inv(params.σ)))))

    # change status and set infective responsible
    S.status = :E
    S.infected_by = I

    # calculate future trajectory
    tI = tE + rand(samplers.α)
    tR = tI + rand(samplers.γ)
    tS = tR + rand(samplers.σ)

    S.t_enter_E = tE
    S.t_enter_I = tI
    S.t_enter_R = tR
    S.t_enter_S = tS

    nothing
end #fn


"""
    make_contacts!(WCW, X, t, dt, params)

Make contacts between individuals in `X` happen, note contacts in `WCW`.
"""
function make_contacts!(WCW, X, t, dt, params)
	N = length(X)
	contact_rate = params.contact_rate
	p_transmission = params.β / contact_rate

	sampler_person = sampler(DiscreteUniform(1,N))
	sampler_num = sampler(Poisson(contact_rate * dt))
	sampler_bern = sampler(Bernoulli(p_transmission))
	samplers = (α = sampler(Exponential(inv(params.α))),
				γ = sampler(Exponential(inv(params.γ))),
    			σ = sampler(Exponential(inv(params.σ))))
	for x in X
		x.status != :I && continue

		num_contacts = rand(sampler_num)
		for i in 1:num_contacts
			# sample with replacement
			# there's no reason why you can't spend multiple contacts on the same
			# individual and increase probability of transmission
			c = rand(sampler_person)

			# note contact
			WCW[c, i] = params.track_steps

			if X[c].status == :S && rand(sampler_bern)
				expose!(X[c], i, t, params, samplers)
			end #if
		end #for
	end # for

	nothing
end #fn


"""
    update_statuses!(X, t)

Update statuses for all individuals in `X` at time `t`.
"""
function update_statuses!(X, t)
    for x in X
        x.status == :S && continue

		if x.status == :E && t ≥ x.t_enter_I
			x.status = :I
		end #if

		if x.status == :I && t ≥ x.t_enter_R
			x.status = :R
		end #if

		if x.status == :R && t ≥ x.t_enter_S
			x.status = :S
		end #if

    end #for

    nothing
end #fn



"""
    get_var(X, v)

Make a list of indices of all individuals in `X` with status `v`.
"""
get_var(X, v) = findall(x -> x.status == v, X)

# get_susceptibles(X) = findall(x -> x.status == :S, X)
# get_infectives(X) = findall(x -> x.status == :I, X)


"""
    get_sum(X, v)

Get the total no. of individuals in `X` with status `v`.
"""
function get_sum(X, v)
	total = 0
	for x in X
		total += x.status == v
	end #for
	total
end #fn
