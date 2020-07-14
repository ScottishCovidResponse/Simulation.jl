using JLSO
using SparseArrays

"""
    AbstractEpiSystem

Abstract supertype for all disease system types.
"""
abstract type AbstractEpiSystem{Part <: AbstractEpiEnv, EL <: EpiList, TR <: AbstractTraitRelationship} <: AbstractMetacommunity{Float64, Matrix{Int64},
                                    Matrix{Float64}, EL, Part} end


mutable struct EpiCache
  virusdecay::Array{Float64, 2}
  virusmigration::Array{Float64, 2}
  valid::Bool
end

mutable struct EpiLookup
  homelookup::SparseMatrixCSC{Float64,Int64}
  worklookup::SparseMatrixCSC{Float64,Int64}
end

"""
    EpiSystem{EE <: AbstractEpiEnv, EL <: EpiList, ER <: AbstractRelationship} <: AbstractEpiSystem{EE, EL, ER}

EpiSystem houses information on different disease classes, `epilist`, the environment, `epienv`, and their relationship to one another, `relationship`.

See `help?>plot_epidynamics` and `help?>plot_epiheatmaps` for relevant plotting functions.
"""
mutable struct EpiSystem{U <: Integer, EE <: AbstractEpiEnv, EL <: EpiList, ER <: AbstractTraitRelationship} <: AbstractEpiSystem{EE, EL, ER}
  abundances::EpiLandscape{U}
  epilist::EL
  epienv::EE
  ordinariness::Union{Matrix{Float64}, Missing}
  relationship::ER
  lookup::EpiLookup
  cache::EpiCache

  function EpiSystem{U, EE, EL, ER}(abundances::EpiLandscape{U},
    epilist::EL, epienv::EE, ordinariness::Union{Matrix{Float64}, Missing}, relationship::ER, lookup::EpiLookup, cache::EpiCache) where {U <: Integer, EE <:
     AbstractEpiEnv,
    EL <: EpiList, ER <: AbstractTraitRelationship}
    new{U, EE, EL, ER}(abundances, epilist, epienv, ordinariness, relationship, lookup, cache)
  end
end

function EpiSystem(popfun::F, epilist::EpiList, epienv::GridEpiEnv,
    rel::AbstractTraitRelationship, intnum::U) where {F<:Function, U <: Integer}

  # Create matrix landscape of zero abundances
  ml = emptyepilandscape(epienv, epilist, intnum)
  # Populate this matrix with species abundances
  popfun(ml, epilist, epienv, rel)
  # Create lookup table of all moves and their probabilities
  home_lookup = genlookups(epienv, epilist.human.movement.home)
  work_lookup = genlookups(epienv, epilist.human.movement.work)
  lookup = EpiLookup(home_lookup, work_lookup)
  nm = zeros(Float64, size(ml.matrix))
  vm = zeros(Float64, size(ml.matrix))
  EpiSystem{U, typeof(epienv), typeof(epilist), typeof(rel)}(ml, epilist, epienv, missing, rel, lookup, EpiCache(nm, vm, false))
end

function EpiSystem(epilist::EpiList, epienv::GridEpiEnv, rel::AbstractTraitRelationship, intnum::U = Int64(1)) where U <: Integer
    return EpiSystem(populate!, epilist, epienv, rel, intnum)
end

function EpiSystem(epilist::EpiList, epienv::GridEpiEnv, rel::AbstractTraitRelationship, initial_population, intnum::U = Int64(1)) where U <: Integer
    if size(initial_population) != size(epienv.active)
        msg = "size(initial_population)==$(size(initial_population)) != " *
            "size(epienv.active)==$(size(epienv.active))"
        throw(DimensionMismatch(msg))
    end
    epienv.active .&= .!_inactive.(initial_population)
    epi = EpiSystem(epilist, epienv, rel, intnum)
    # Add in the initial susceptible population
    idx = findfirst(epilist.human.names .== "Susceptible")
    if idx == nothing
        msg = "epilist has no Susceptible category. epilist.names = $(epilist.human.names)"
        throw(ArgumentError(msg))
    end
    # Modify active cells based on new population
    initial_population = convert_population(initial_population, intnum)
    epi.abundances.grid[idx, :, :] .+= initial_population
    return epi
end

"""
    isapprox(epi_1::AbstractEpiSystem, epi_2::AbstractEpiSystem; kwargs...)

Compare two `EpiSystem`s for approximate equality. Specifically, compares the
`EpiLandscape`s of the two systems.

## Keyword arguments
- Anything to pass to `Base.isapprox`.

!!! note
    You may want to pass in `atol` or `rtol` to loosen the equality tolerance.
"""
function Base.isapprox(epi_1::AbstractEpiSystem, epi_2::AbstractEpiSystem; kwargs...)
    return isapprox(epi_1.abundances, epi_2.abundances; kwargs...)
end

save(path::String, system::EpiSystem) = JLSO.save(path, :episystem => system)
load(path::String, obj_type::Type{EpiSystem}) = JLSO.load(path)[:episystem]

function getsize(epi::AbstractEpiSystem)
  return _getsize(epi.epienv.habitat)
end

function getgridsize(epi::AbstractEpiSystem)
  return _getgridsize(epi.epienv.habitat)
end

function getdimension(epi::AbstractEpiSystem)
    return _getdimension(epi.epienv.habitat)
end


function gettraitrel(epi::AbstractEpiSystem)
  return epi.relationship
end

function gethabitat(epi::AbstractEpiSystem)
  return epi.epienv.habitat
end

import Diversity.API: _getabundance
function _getabundance(epi::AbstractEpiSystem, input::Bool)
    if input
        return epi.abundances.matrix
    else
        return _calcabundance(_gettypes(epi), epi.abundances.matrix / sum(epi.abundances.matrix))[1]
    end
end

import Diversity.API: _getmetaabundance
function _getmetaabundance(epi::AbstractEpiSystem)
  return sumoversubcommunities(epi, _getabundance(epi))
end


import Diversity.API: _getpartition
function _getpartition(epi::AbstractEpiSystem)
  return epi.epienv
end
import Diversity.API: _gettypes
function _gettypes(epi::AbstractEpiSystem)
    return epi.epilist
end
import Diversity.API: _getordinariness!
function _getordinariness!(epi::AbstractEpiSystem)
    if ismissing(epi.ordinariness)
        relab = getabundance(epi, false)
        epi.ordinariness = _calcordinariness(epi.epilist, relab)
    end
    return epi.ordinariness
end

import Diversity.API._getscale
function _getscale(epi::AbstractEpiSystem)
    return _calcabundance(_gettypes(epi), getabundance(epi, false))[2]
end

function invalidatecaches!(epi::AbstractEpiSystem)
    epi.ordinariness = missing
    epi.cache.virusdecay .= 0
    epi.cache.virusmigration .= 0
    epi.cache.valid = false
end

function getdispersaldist(epi::AbstractEpiSystem, sp::Int64)
  dist = epi.epilist.human.movement.home.kernels[sp].dist
  return dist
end
function getdispersaldist(epi::AbstractEpiSystem, sp::String)
  num = Compat.findall(epi.epilist.human.names.==sp)[1]
  getdispersaldist(epi, num)
end

function getdispersalvar(epi::AbstractEpiSystem, sp::Int64)
    var = (epi.epilist.human.movement.home.kernels[sp].dist)^2 * pi / 4
    return var
end
function getdispersalvar(epi::AbstractEpiSystem, sp::String)
    num = Compat.findall(epi.epilist.human.names.==sp)[1]
    getdispersalvar(epi, num)
end

function getlookup(epi::AbstractEpiSystem, id::Int64, movetype::String)
    if movetype == "home"
        return epi.lookup.homelookup[id, :]
    elseif movetype == "work"
        return epi.lookup.worklookup[id, :]
    else
        return error("No other movement types currently implemented")
    end
end

function getlookup(epi::AbstractEpiSystem, id::Int64)
    return epi.lookup.homelookup[id, :], epi.lookup.worklookup[id, :]
end

function genlookups(epienv::AbstractEpiEnv, mov::Commuting)
    total_size = (size(epienv.active, 1) * size(epienv.active, 2))
    Is = Int64.(mov.home_to_work[!, :from])
    Js = Int64.(mov.home_to_work[!, :to])
    Vs = mov.home_to_work[!, :count]
    work = sparse(Is, Js, Vs, total_size, total_size)
    dropzeros!(work)
    for i in work.rowval
        work[i, :] ./= sum(work[i, :])
    end
    return sparse(Is, Js, Vs, total_size, total_size)
end
function genlookups(epienv::GridEpiEnv, mov::AlwaysMovement)
    total_size = (size(epienv.active, 1) * size(epienv.active, 2))
    grid_locs = 1:total_size
    activity = epienv.active[1:end]
    grid_locs = grid_locs[activity]
    xys = convert_coords.(grid_locs, size(epienv.active, 2))
    grid_size = _getgridsize(epienv.habitat)
    sd = [(2 .* k.dist) ./ sqrt(pi) for k in mov.kernels][activity]
    relsize =  grid_size ./ sd
    thresh = [k.thresh for k in mov.kernels][activity]
    grid_size /= unit(grid_size)
    res = map((i, r, t) -> Simulation.genlookups(i, grid_locs, xys, grid_size, r, t, epienv), grid_locs, relsize, thresh)
    Is = vcat([fill(grid_locs[r], length(res[r][1])) for r in eachindex(res)]...)
    Js = vcat([r[1] for r in res]...)
    Vs = vcat([r[2] for r in res]...)
    return sparse(Is, Js, Vs, total_size, total_size)
end

function genlookups(from::Int64, to::Vector{Int64}, xys::Array{Tuple{Int64,Int64},1}, grid_size::Float64, relsize::Float64, thresh::Float64, epienv::GridEpiEnv)
    x, y = xys[to .== from][1]
    maxX = ceil(Int64, x + grid_size * 2/relsize); minX = ceil(Int64, x - grid_size * 2/relsize)
    maxY = floor(Int64, y + grid_size * 2/relsize); minY = floor(Int64, y - grid_size * 2/relsize)
    keep = [(i[1] <= maxX) & (i[2] <= maxY) & (i[1] >= minX) & (i[2] >= minY) for i in xys]
    to = to[keep]
    probs = [_lookup((x = x, y = y), (x = i[1], y = i[2]), relsize, _gaussian_disperse) for i in xys[keep]]
    keep = probs .> thresh
    probs = probs[keep]
    probs ./= sum(probs)
    return to[keep], probs
end

function _lookup(from::NamedTuple, to::NamedTuple, relSquareSize::Float64, dispersalfn::F) where {F<:Function}
    return calc_prob = hcubature(r -> dispersalfn(r),
      [from.x - relSquareSize, from.y - relSquareSize, to.x - relSquareSize, to.y - relSquareSize],
      [from.x, from.y, to.x, to.y],
      maxevals= 100, rtol = 0.01)[1] / relSquareSize^2
end
