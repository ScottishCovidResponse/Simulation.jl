module Persons

export Person

mutable struct Person
	status::Symbol
	t_enter_E::Float64
	t_enter_I::Float64
	t_enter_R::Float64
	t_enter_S::Float64
	t_test::Float64
	t_isolated::Float64
	infected_by::Int
	alerted::Float64
	alerted_by::Int
end #struct

function Person(;
	status = :S,
	t_enter_E = Inf,
	t_enter_I = Inf,
	t_enter_R = Inf,
	t_enter_S = Inf,
	t_test = Inf,
	t_isolated = Inf,
	infected_by = -1,
	alerted = Inf,
	alerted_by = -1)

	Person(status, t_enter_E, t_enter_I, t_enter_R, t_enter_S, t_test, t_isolated, infected_by, alerted, alerted_by)
end #fn

end #module
