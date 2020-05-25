
"""
    EpiList{TR <: AbstractTraits,
                 MO <: AbstractMovement,
                 T <: AbstractTypes,
                 P <: AbstractParams} <: AbstractTypes
Epi list houses all disease class specific information including trait information and
movement types.
"""
mutable struct EpiList{
    TR <: AbstractTraits,
    MO <: AbstractMovement,
    T <: AbstractTypes,
    P <: AbstractParams,
} <: AbstractTypes

    disease_classes::Vector{String}
    risk_factors::NamedTuple
    traits::TR
    abun::Vector{Int64}
    types::T
    movement::MO
    params::P

    function EpiList{TR, MO, T, P}(
        disease_classes::Vector{String},
        risk_factors::NamedTuple,
        traits::TR,
        abun::Vector{Int64},
        types::T,
        movement::MO,
        params::P
    ) where {TR<:AbstractTraits, MO<:AbstractMovement, T<:AbstractTypes, P<:AbstractParams}
        return new{TR, MO, T, P}(
            disease_classes, risk_factors, traits, abun, types, movement, params,
        )
    end

    function EpiList{TR, MO, T, P}(
        traits::TR, abun::Vector{Int64}, types::T, movement::MO, params::P
    ) where {
        TR <: AbstractTraits,
        MO <: AbstractMovement,
        T <: AbstractTypes,
        P <: AbstractParams
    }
        disease_classes = map(x -> "$x", 1:length(abun))
        risk_factors = NamedTuple() # Assume none
        new{TR, MO, T, P}(
            disease_classes, risk_factors, traits, abun, types, movement, params
        )
    end
end

import Diversity.API: _gettypenames
function _gettypenames(el::EpiList, input::Bool)
    return _gettypenames(el.types, input)
end
import Diversity.API: _counttypes
function _counttypes(el::EpiList, input::Bool)
    return _counttypes(el.types, input)
end
import Diversity.API: _calcsimilarity
function _calcsimilarity(el::EpiList, a::AbstractArray)
    return _calcsimilarity(el.types, a)
end
import Diversity.API: floattypes
function floattypes(::EpiList)
    return Set([Float64])
end
import Diversity.API: _calcordinariness
function _calcordinariness(el::EpiList, a::AbstractArray)
    _calcordinariness(el.types, a, one(eltype(a)))
end
import Diversity.API: _calcabundance
function _calcabundance(el::EpiList, a::AbstractArray)
  return _calcabundance(el.types, a)
end
import Diversity.API._getdiversityname
function _getdiversityname(el::EpiList)
    return _getdiversityname(el.types)
end

function SIS(
    traits::TR,
    abun::Vector{Int64},
    movement::MO,
    params::P,
    risk_factors = NamedTuple(),
) where {TR <: AbstractTraits, MO <: AbstractMovement, P <: AbstractParams}

    disease_classes = ["Virus", "Susceptible", "Infected", "Dead"]
    total_classes = if length(risk_factors) == 0
        length(disease_classes)
    else
        length(disease_classes) * sum(values(a))
    end
    types = UniqueTypes(length(total_classes))
    length(abun) == length(total_classes) || throw(
        DimensionMismatch("Abundance vector doesn't match number of disease classes")
    )

    return EpiList{typeof(traits), typeof(movement), typeof(types), typeof(params)}(
        disease_classes, risk_factors, traits, abun, types, movement, params
    )
end

function SIR(traits::TR, abun::Vector{Int64},
    movement::MO, params::P, age_categories::Int64 = 1) where {TR <: AbstractTraits,
        MO <: AbstractMovement, P <: AbstractParams}

    names = ["Susceptible", "Infected", "Recovered", "Dead"]
    new_names = [ifelse(i == 1, "$j", "$j$i") for i in 1:age_categories, j in names]
    new_names =  ["Virus"; new_names[1:end]]
    types = UniqueTypes(length(new_names))
    length(abun) == length(new_names) || throw(DimensionMismatch("Abundance vector doesn't match number of disease classes"))
  EpiList{typeof(traits), typeof(movement), typeof(types), typeof(params)}(new_names, traits, abun, types, movement, params)
end

function SEIR(traits::TR, abun::Vector{Int64},
    movement::MO, params::P, age_categories::Int64 = 1) where {TR <: AbstractTraits,
        MO <: AbstractMovement, P <: AbstractParams}

    names = ["Susceptible", "Exposed", "Infected", "Recovered", "Dead"]
    new_names = [ifelse(i == 1, "$j", "$j$i") for i in 1:age_categories, j in names]
    new_names =  ["Virus"; new_names[1:end]]
    types = UniqueTypes(length(new_names))
    length(abun) == length(new_names) || throw(DimensionMismatch("Abundance vector doesn't match number of disease classes"))
  EpiList{typeof(traits), typeof(movement), typeof(types), typeof(params)}(new_names, traits, abun, types, movement, params)
end

function SEIRS(traits::TR, abun::Vector{Int64},
    movement::MO, params::P, age_categories::Int64 = 1) where {TR <: AbstractTraits,
        MO <: AbstractMovement, P <: AbstractParams}

    return SEIR(traits, abun, movement, params, age_categories)
end

function SEI2HRD(traits::TR, abun::Vector{Int64},
    movement::MO, params::P, age_categories::Int64 = 1) where {TR <: AbstractTraits,
        MO <: AbstractMovement, P <: AbstractParams}

    names = ["Susceptible", "Exposed", "AsymptomaticInfected", "SymptomaticInfected", "Hospitalised", "Recovered", "Dead"]
    new_names = [ifelse(i == 1, "$j", "$j$i") for i in 1:age_categories, j in names]
    new_names =  ["Virus"; new_names[1:end]]
    types = UniqueTypes(length(new_names))
    length(abun) == length(new_names) || throw(DimensionMismatch("Abundance vector doesn't match number of disease classes"))
    size(params.transition, 1) == (length(new_names) - 1) || throw(DimensionMismatch("Transition matrix does not have the correct number of classes"))
  EpiList{typeof(traits), typeof(movement), typeof(types), typeof(params)}(new_names, traits, abun, types, movement, params)
end
