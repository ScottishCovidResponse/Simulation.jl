import Base.eltype

"""
    AbstractRequirement{Energy}

Abstract supertype for all species energy requirement types, parameterised by
the type(s) of energy required `Energy`.
"""
abstract type AbstractRequirement{Energy} end

function eltype{Energy}(::AbstractRequirement{Energy})
  return Energy
end

"""
    SimpleRequirement <: AbstractRequirement{Float64}

A simple energy requirement is a single float for each species.
"""
mutable struct SimpleRequirement <: AbstractRequirement{Float64}
  energy::Vector{Float64}
end

"""
    AbstractBudget

Abstract supertype for all budget types
"""
abstract type AbstractBudget{Requirement} end

function eltype{Energy}(::AbstractBudget{Energy})
  return Energy
end

function countsubcommunities(ab::AbstractBudget)
  return _countsubcommunities(ab)
end

"""
    SimpleBudget <: AbstractBudget{Float64}

This budget type has a matrix of floats, representing the energy budget of each
subcommunity in the abiotic environment.
"""
mutable struct SimpleBudget <: AbstractBudget{Float64}
  matrix::Matrix{Float64}
end

function _countsubcommunities(bud::SimpleBudget)
  return length(bud.matrix)
end
