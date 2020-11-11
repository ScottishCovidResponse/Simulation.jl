using DataFrames

mutable struct EpiTransition
    virus_effects::Array{Float64, 2}
    transition_effects::Array{Float64, 2}
    transitions_from::Vector{Int64}
    transitions_to::Vector{Int64}
end

function EpiTransition(epilist::EpiList, transitions_from::Vector{Int64}, transitions_to::Vector{Int64})
    ve = fill(1.0, size(epilist.params.transition_force))
    te = fill(1.0, size(epilist.params.transition))
    return EpiTransition(ve, te, transitions_from, transitions_to)
end

function EpiTransition(epilist::EpiList)
    ve = fill(1.0, size(epilist.params.transition_force))
    te = fill(1.0, size(epilist.params.transition))
    return EpiTransition(ve, te, [1], [1])
end
