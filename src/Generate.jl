using StatsBase
using Query
using Compat

"""
    get_neighbours(mat::Matrix, x_coord::Int64, y_coord::Int64, chess::Int64=4)

Function to get the neighbours of a grid square in a matrix in 4 or 8 directions
"""
function get_neighbours(mat::Matrix, x_coord::Int64, y_coord::Int64, chess::Int64=4)
  # Calculate dimensions
  dims=size(mat)
  x_coord <= dims[1] && y_coord <= dims[2] || error("Coordinates outside grid")
  # Include 4 directions
  if chess==4
    neighbour_vec=[x_coord y_coord-1; x_coord y_coord+1; x_coord-1 y_coord;
     x_coord+1 y_coord]
  # Include 8 directions
  elseif chess==8
    neighbour_vec=[x_coord y_coord-1; x_coord y_coord+1; x_coord-1 y_coord;
     x_coord+1 y_coord; x_coord-1 y_coord-1; x_coord-1 y_coord+1;
      x_coord+1 y_coord-1; x_coord+1 y_coord+1]
  else
    # Give error if other number chosen than 4 or 8
    error("Can only calculate neighbours in 4 or 8 directions")
  end
  # Remove answers outside of the dimensions of the matrix
  remove=vcat(mapslices(all, [neighbour_vec.>=1 neighbour_vec[:,1].<=
    dims[1] neighbour_vec[:,2].<=dims[2]], dims=2)...)
  neighbour_vec=neighbour_vec[remove,:]
  neighbour_vec
end
function get_neighbours(mat::Matrix, x_coord::Array{Int64,1},
     y_coord::Array{Int64,1}, chess::Int64=4)
     neighbours  =map(n -> get_neighbours(mat, x_coord[n], y_coord[n], chess),
      eachindex(x_coord))
      return vcat(neighbours...)
end
"""
    update!(eco::Ecosystem,  birth::Float64, death::Float64,
       l::Float64, s::Float64, timestep::Real)
Function to update a Ecosystem after one timestep. It takes in parameters of
birth, death rates and longevity of species (l & s) to generate the abundances
of the species stochastically. Movement takes place across the landscape via
movement rates defined in the ecosystem.
"""
function update!(eco::Ecosystem, timestep::Unitful.Time)
    # Calculate dimenions of habitat and number of species
    dims = _countsubcommunities(eco.abenv.habitat)
    spp = size(eco.abundances.grid,1)
    net_migration = eco.cache.netmigration
    params = eco.spplist.params
      # Loop through grid squares
      for i in 1:dims
          # Get the overall energy budget of that square
          width = getdimension(eco)[1]
          (x, y) = convert_coords(i, width)
          if (eco.abenv.active[x, y] && sum(eco.abundances.matrix[:, i])!=0)
             birth_energy, death_energy = energy(eco, eco.abenv.budget, i,
              spp)
              currentabun = eco.abundances.matrix[:, i]
              # Loop through species in chosen square
              for j in 1:spp

                # Calculate effective rates
                birthprob = params.birth[j] * timestep * birth_energy[j]
                deathprob = params.death[j] * timestep * death_energy[j]

                # Put probabilities into 0 - 1
                newbirthprob = 1.0 - exp(-birthprob)
                newdeathprob = 1.0 - exp(-deathprob)

                # Calculate how many births and deaths
                births = jbinom(1, currentabun[j], newbirthprob)[1]
                deaths = jbinom(1, currentabun[j], newdeathprob)[1]
                # Update population
                eco.abundances.matrix[j, i] += (births - deaths)

                # Perform gaussian movement
                move!(eco, eco.spplist.movement, i, j, net_migration, births)
            end
        end
    end
    eco.abundances.matrix .= eco.abundances.matrix .+ net_migration
    eco.cache.netmigration .= 0
    resetcache!(eco)
    # Update environment
    habitatupdate!(eco, timestep)
    budgetupdate!(eco, timestep)
end

function energy(eco::Ecosystem, bud::AbstractBudget, i::Int64, spp::Int64)
    if typeof(eco.spplist.params) <: NoGrowth
        birth_energy = Vector{Float64}(spp); death_energy = Vector{Float64}(spp)
        fill!(birth_energy, 0.0); fill!(death_energy, 0.0)
    else
        width = getdimension(eco)[1]
        (x, y) = convert_coords(i, width)
        params = eco.spplist.params
        K = ustrip.(getbudget(eco)[x, y])
        # Get abundances of square we are interested in
        currentabun = eco.abundances.matrix[:, i]

        # Get energy budgets of species in square
        ϵ̄ = ustrip.(eco.spplist.requirement.energy)
        E = sum(convert(Vector{Float64}, currentabun) .* ϵ̄)
        # Traits
        ϵ̄real = copy(ϵ̄)
        birth_energy = Vector{Float64}(Compat.undef, spp); death_energy = Vector{Float64}(Compat.undef, spp)
        for k in 1:spp
          ϵ̄real[k] = ϵ̄[k]/traitfun(eco, i, k)
          # Alter rates by energy available in current pop & own requirements
          birth_energy[k] = ϵ̄[k]^-params.l * ϵ̄real[k]^-params.s * min(K/E, params.boost)
          death_energy[k] = ϵ̄[k]^-params.l * ϵ̄real[k]^params.s * (E / K)
        end
    end
    return birth_energy, death_energy
end
function energy(eco::Ecosystem, bud::BudgetCollection2, i::Int64, spp::Int64)
     width = getdimension(eco)[1]
     (x, y) = convert_coords(i, width)
     params = eco.spplist.params
    K1 = ustrip.(_getbudget(eco.abenv.budget, :b1)[x, y])
    K2 = ustrip.(_getbudget(eco.abenv.budget, :b2)[x, y])
    # Get abundances of square we are interested in
    currentabun = eco.abundances.matrix[:, i]

    # Get energy budgets of species in square
    ϵ̄1 = ustrip.(eco.spplist.requirement.r1.energy)
    E1 = sum(convert(Vector{Float64}, currentabun) .* ϵ̄1)
    ϵ̄2 = ustrip.(eco.spplist.requirement.r2.energy)
    E2 = sum(convert(Vector{Float64}, currentabun) .* ϵ̄2)
    # Traits
    ϵ̄real1 = copy(ϵ̄1)
    ϵ̄real2 = copy(ϵ̄2)
    birth_energy = Vector{Float64}(Compat.undef, spp); death_energy = Vector{Float64}(Compat.undef, spp)
    for k in 1:spp
      ϵ̄real1[k] = ϵ̄1[k]/traitfun(eco, i, k)
      ϵ̄real2[k] = ϵ̄2[k]/traitfun(eco, i, k)
      # Alter rates by energy available in current pop & own requirements
      birth_energy[k] = (ϵ̄1[k]^-params.l + ϵ̄2[k]^-params.l) * (ϵ̄real1[k]^-params.s +
         ϵ̄real2[k]^-params.s) * min(K1/E1, K2/E2, params.boost)
      death_energy[k] = (ϵ̄1[k]^-params.l + ϵ̄2[k]^-params.l) * (ϵ̄real1[k]^params.s +
        ϵ̄real2[k]^params.s) * max(E1 / K1, E2/K2)
    end
    return birth_energy, death_energy
end


function convert_coords(i::Int64, width::Int64)
  x = ((i - 1) % width) + 1
  y = div((i - 1), width)  + 1
  return (x, y)
end
function convert_coords(i::Array{Int64, 1}, width::Int64)
  x = ((i .- 1) .% width) .+ 1
  y = div.((i .- 1), width)  .+ 1
  return (x, y)
end
function convert_coords(x::Int64, y::Int64, width::Int64)
  i = x + width * (y - 1)
  return i
end

function convert_coords(x::Array{Int64, 1}, y::Array{Int64, 1}, width::Int64)
  i = x .+ (width .* (y .- 1))
  return i
end
function calc_lookup_moves(bound::NoBoundary, x::Int64, y::Int64, spp::Int64, eco::Ecosystem, abun::Int64)
  lookup = eco.lookup[spp]
  maxX = getdimension(eco)[1] - x
  maxY = getdimension(eco)[2] - y
  # Can't go over maximum dimension
  valid = (lookup.x .> -x) .& (lookup.y .> -y) .&
   (lookup.x .<= maxX) .& (lookup.y .<= maxY)
  for i in eachindex(valid)
      if valid[i]
          valid[i] = valid[i] & (eco.abenv.active[lookup.x[i] .+ x,
          lookup.y[i] .+ y])
      end
  end
  lookup.pnew[.!valid] .= 0.0
  lookup.pnew[valid] = lookup.p[valid]
  lookup.pnew ./= sum(lookup.pnew)
  multinom = Multinomial(abun, lookup.pnew)
  draw = rand(multinom)
  lookup.moves .= draw
end

function calc_lookup_moves(bound::Cylinder, x::Int64, y::Int64, spp::Int64, eco::Ecosystem, abun::Int64)
  lookup = eco.lookup[spp]
  maxX = getdimension(eco)[1] - x
  maxY = getdimension(eco)[2] - y
  lookup.x[lookup.x .<= -x] .+= Simulation.getdimension(eco)[1]
  lookup.x[lookup.x .> maxX] .-= Simulation.getdimension(eco)[1]

  valid = (lookup.y .> -y) .& (lookup.y .<= maxY)
  for i in eachindex(lookup.x)
      if valid[i]
          valid[i] = valid[i] & (eco.abenv.active[lookup.x[i] .+ x,
            lookup.y[i] .+ y])
      end
  end

  lookup.pnew[.!valid] .= 0.0
  lookup.pnew[valid] = lookup.p[valid]
  lookup.pnew ./= sum(lookup.pnew)
  multinom = Multinomial(abun, lookup.pnew)
  draw = rand(multinom)
  lookup.moves .= draw
end

function calc_lookup_moves(bound::Torus, x::Int64, y::Int64, spp::Int64, eco::Ecosystem, abun::Int64)
  lookup = eco.lookup[spp]
  maxX = Simulation.getdimension(eco)[1] - x
  maxY = Simulation.getdimension(eco)[2] - y
  # Can't go over maximum dimension
  lookup.x[lookup.x .<= -x] += Simulation.getdimension(eco)[1]
  lookup.y[lookup.y .<= -y] += Simulation.getdimension(eco)[2]
  lookup.x[lookup.x .> maxX] -= Simulation.getdimension(eco)[1]
  lookup.y[lookup.y .> maxX] -= Simulation.getdimension(eco)[2]

  valid = (lookup.x .> -x) .& (lookup.y .> -y) .&
   (lookup.x .<= maxX) .& (lookup.y .<= maxY)
  for i in eachindex(valid)
      if valid[i]
          valid[i] = valid[i] & (eco.abenv.active[lookup.x[i] .+ x,
          lookup.y[i] .+ y])
      end
  end

  lookup.pnew[.!valid] .= 0.0
  lookup.pnew[valid] = lookup.p[valid]
  lookup.pnew ./= sum(lookup.pnew)
  multinom = Multinomial(abun, lookup.pnew)
  draw = rand(multinom)
  lookup.moves .= draw
end
"""
    move!(i::Int64, spp::Int64, eco::Ecosystem, grd::Array{Int64, 2})

Function to calculate the movement of species `spp` from a given position in the
landscape `i`, using the lookup table found in the Ecosystem and updating the
movement patterns on a grid, `grd`. Optionally, a number of births can be
provided, so that movement only takes place as part of the birth process, instead
of the entire population
"""
function move!(eco::Ecosystem, ::AlwaysMovement, i::Int64, spp::Int64,
  grd::Array{Int64, 2}, ::Int64)

  width = getdimension(eco)[1]
  (x, y) = convert_coords(i, width)
  full_abun = eco.abundances.matrix[spp, i]
  calc_lookup_moves(getboundary(eco.spplist.movement), x, y, spp, eco, full_abun)
  # Lose moves from current grid square
  grd[spp, i] = grd[spp, i] - sum(eco.lookup[spp].moves)
  valid = eco.lookup[spp].pnew .> 0.0
  # Map moves to location in grid
  mov = eco.lookup[spp].moves[valid]
  locs = convert_coords((eco.lookup[spp].x .+ x), (eco.lookup[spp].y .+ y), width)[valid]
  for i in eachindex(locs)
     grd[spp, locs[i]] += mov[i]
  end
  return eco
end
"""
    move!(i::Int64, spp::Int64, eco::Ecosystem, grd::Array{Int64, 2})

Function to calculate the movement of species `spp` from a given position in the
landscape `i`, using the lookup table found in the Ecosystem and updating the
movement patterns on a grid, `grd`. Optionally, a number of births can be
provided, so that movement only takes place as part of the birth process, instead
of the entire population
"""
function move!(eco::Ecosystem, ::NoMovement, i::Int64, spp::Int64,
  grd::Array{Int64, 2}, ::Int64)
  return eco
end


    function move!(eco::Ecosystem, ::BirthOnlyMovement, i::Int64, spp::Int64,
        grd::Array{Int64, 2}, births::Int64)
      width = getdimension(eco)[1]
      (x, y) = convert_coords(i, width)
      table = calc_lookup_moves(getboundary(eco.spplist.movement), x, y, spp, eco, births)
      # Lose moves from current grid square
      grd[spp, i] = grd[spp, i] - sum(eco.lookup[spp].moves)
      valid = eco.lookup[spp].pnew .> 0.0
      # Map moves to location in grid
      mov = eco.lookup[spp].moves[valid]
      locs = convert_coords((eco.lookup[spp].x .+ x), (eco.lookup[spp].y .+ y), width)[valid]
      for i in eachindex(locs)
         grd[spp, locs[i]] += mov[i]
      end
      return eco
    end



"""
    populate!(ml::GridLandscape, spplist::SpeciesList,
                   abenv::AbstractAbiotic, traits::Bool)

Function to populate a grid landscape given the abundances found in species list
and whether or not to include traits.
"""
function populate!(ml::GridLandscape, spplist::SpeciesList,
                   abenv::AbstractAbiotic)
  # Calculate size of habitat
  dim = _getdimension(abenv.habitat)
  len = dim[1] * dim[2]
  grid = collect(1:len)
  # Set up copy of budget
      b = reshape(copy(_getbudget(abenv.budget)), size(grid))
      units = unit(b[1])
      activity = reshape(copy(abenv.active), size(grid))
      b[.!activity] .= 0.0 * units
      # Loop through species
      for i in eachindex(spplist.abun)
        # Get abundance of species
        abun = spplist.abun[i]
        # Loop through individuals
          while abun>0
            # Randomly choose position on grid (weighted)
            pos = sample(grid[b .> (0 * units)])
          # Add individual to this location
          ml.matrix[i, pos] = ml.matrix[i, pos] .+ 1
          abun = abun .- 1
          b[pos] = b[pos] .- spplist.requirement.energy[i]
        end
      end
end
function populate!(ml::GridLandscape, spplist::SpeciesList,
                   abenv::GridAbioticEnv{H, BudgetCollection2{B1, B2}} where {H <:
                    AbstractHabitat, B1 <: AbstractBudget, B2 <: AbstractBudget})
    # Calculate size of habitat
    dim = _getdimension(abenv.habitat)
    len = dim[1] * dim[2]
    grid = collect(1:len)
    # Set up copy of budget
    b1 = reshape(copy(_getbudget(abenv.budget, :b1)), size(grid))
    b2 = reshape(copy(_getbudget(abenv.budget, :b2)), size(grid))
    units1 = unit(b1[1])
    units2 = unit(b2[1])
    activity = reshape(copy(abenv.active), size(grid))
    b1[.!activity] = 0.0 * units1
    b2[.!activity] = 0.0 * units2
    # Loop through species
    for i in eachindex(spplist.abun)
        # Get abundance of species
        abun = spplist.abun[i]
        # Loop through individuals
        while abun>0
            # Randomly choose position on grid (weighted)
            pos = sample(grid[(b1 .> (0 * units1)) .& (b2 .> (0 * units2))])
            # Add individual to this location
            ml.matrix[i, pos] = ml.matrix[i, pos] .+ 1
            abun = abun .- 1
            b1[pos] = b1[pos] .- spplist.requirement.r1.energy[i]
            b2[pos] = b2[pos] .- spplist.requirement.r2.energy[i]
        end
    end
end

function simplepopulate!(ml::GridLandscape, spplist::SpeciesList,
                   abenv::AbstractAbiotic)
  # Calculate size of habitat
  dim = _getdimension(abenv.habitat)
  len = dim[1] * dim[2]
  grid = collect(1:len)
  # Set up copy of budget
  b = reshape(copy(_getbudget(abenv.budget)), size(grid))
  units = unit(b[1])
  activity = reshape(copy(abenv.active), size(grid))
  b[.!activity] .= 0.0 * units
  # Loop through species
  for i in eachindex(spplist.abun)
      if spplist.native[i]
        # Get abundance of species
        abun = spplist.abun[i]
        pos = sample(grid[b .> (0 * units)], abun)
        # Add individual to this location
        ml.matrix[i, pos] = ml.matrix[i, pos] .+ 1
     end
   end
end
"""
    repopulate!(eco::Ecosystem)
Function to repopulate an ecosystem `eco`, with option for including trait
preferences.
"""
function repopulate!(eco::Ecosystem)
  eco.abundances = emptygridlandscape(eco.abenv, eco.spplist)
  eco.spplist.abun = rand(Multinomial(sum(eco.spplist.abun), length(eco.spplist.abun)))
  populate!(eco.abundances, eco.spplist, eco.abenv)
end

"""
    reenergise!(eco::Ecosystem, budget::)
Function to repopulate an ecosystem `eco`, with option for including trait
preferences.
"""
function reenergise!(eco::Ecosystem, budget::Float64, grid::Tuple{Int64, Int64})
    fill!(eco.abenv.budget.matrix, budget/(grid[1]*grid[2]))
end
