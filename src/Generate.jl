using StatsBase

# Function to randomly populate a habitat matrix
function populate(species::Int64, individuals::Int64, habitat::Habitats,
   dist::Distribution= Multinomial(individuals,species))
  # Calculate size of habitat
  dim=size(habitat.matrix)
  # Create empty population matrix of correct dimensions
  P=zeros(Int64,species,dim[1],dim[2])
  # Randomly choose abundances for each species from Multinomial
  abun_vec=rand(dist)
  # Loop through species
  for i in eachindex(abun_vec)
    # Get abundance of species
    abun=abun_vec[i]
      # Loop through individuals
      while abun>0
      # Randomly choose position on grid
      x=rand(1:dim[1])
      y=rand(1:dim[1])
      # Add individual to this location
      P[i,x,y]=P[i,x,y]+1
      abun=abun-1
    end
  end
  # Create MatrixLandscape from P matrix and habitat
  MatrixLandscape(P, habitat)
end

# Function to calculate species richness from an Ecosystem object
function SR(ecosystem::Ecosystem)
  sz=size(ecosystem.partition.abundances,2,3)
  ms=map(x-> sum(ecosystem.partition.abundances[:,x].>0), 1:(sz[1]*sz[2]))#*individuals
  reshape(ms, sz)
end

# Function to create a habitat from a discrete set of types
function create_habitat(dim, types, prop)
  # Weighted sample from the types in the correct dimension
  sample(types, WeightVec(prop), dim)
end

# Function to create a habitat from a discrete set of types according to the
# Saura-Martinez-Millan algorithm (2000)
function percolate!(M::AbstractMatrix, p::Real)
  for i in 1:(length(M))
    if junif(0, 1) < p
      M[i]=1
    end
  end
  # Select clusters and assign types
  count=1
  for x in 1:dim[1]
    for y in 1:dim[2]
      if M[x,y]==1.0

        neighbours=get_neighbours(M, y, x)
        # If no neighbours, sample from clusters randomly
        if isempty(neighbours)
          count=count+1
          M[x,y]=count
        else
        cluster=vcat(mapslices(x->M[x[1],x[2]].==1, neighbours, 2)...)
        already=vcat(mapslices(x->M[x[1],x[2]].>1, neighbours, 2)...)
          if any(already)
            neighbours=neighbours[already,:]
            M[x,y]=M[neighbours[1,1],neighbours[1,2]]
          else
            count=count+1
            neighbours=neighbours[cluster,:]
            M[x,y]=count
            map(i->M[neighbours[i,1],neighbours[i,2]]=count, 1:size(neighbours,1))
          end
        end
      end
    end
  end
  T=Array{String}(dim)
  map(x->T[M.==x]=sample(types, wv), 1:maximum(M))
  # Fill in undefined squares with most frequent neighbour
end

function fill_in!(T, M, wv)
  dim = size(M)
  # Loop through grid of clusters
  for x in 1:dim[1]
    for y in 1:dim[2]
      if M[x,y]==0

        neighbours=get_neighbours(T, y, x, 8)
        # If no neighbours, sample from traits randomly
        if isempty(neighbours)
          T[x,y]=sample(types, wv)
        else
        already=vcat(mapslices(x->isdefined(T,x[1],x[2]), neighbours, 2)...)
        # If any already assigned then sample from most frequent neighbour traits
          if any(already)
            neighbours=neighbours[already,:]

            neighbour_traits=map(i->T[neighbours[i,1],neighbours[i,2]],
             1:size(neighbours,1))
            ind=indmax(map(x->sum(neighbour_traits.==x), types))
            T[x,y]= types[ind]
          # If none are assigned in entire grid already,
          # sample randomly from traits
        elseif all(M.<=1)
            T[x,y]=sample(types, wv)
          # If some are assigned in grid, sample from these
          else
            T[x,y]=sample(T[M.>1])
          end
        end
      end
    end
  end
end
  T

 function random_habitat(dim::Tuple, types, p::Real, A::Vector)
  # Check that the proportion of coverage for each type matches the number
  # of types and that they add up to 1
  length(A)==length(types) || error("There must be an area proportion for each type")
  sum(A)==1 || error("Proportion of habitats must add up to 1")
  # Create weighting from proportion habitats
  wv=weights(A)

  # Create an empty grid of the right dimension
  M=zeros(dim)

  # If the dimensions are too small for the algorithm, just use a weighted sample
  if any(map(x-> dim[x]<=2, 1:2))
    T = create_habitat(dim, types, A)
  else
    # Percolation step
    percolate!(M, p)
    # Select clusters and assign types
    identify_clusters!(M)
    # Create a string grid of the same dimensions
    T=Array{String}(dim)
    # Fill in T with clusters already created
    map(x->T[M.==x]=sample(types, wv), 1:maximum(M))
    # Fill in undefined squares with most frequent neighbour
    fill_in!(T, M, wv)
    T
  end
end

# Function to populate a Niche habitat
function populate(species::Int64, individuals::Int64, habitat::Niches,
  traits::Vector, budget::Budget,
  dist::Distribution= Multinomial(individuals,species))
  # Calculate size of habitat
  dim=size(habitat.matrix)
  grid=collect(1:dim[1]*dim[2])
  # Create empty population matrix of correct dimensions
  P=zeros(Float64,species,dim[1],dim[2])
  # Randomly choose abundances for each species from Multinomial
  abun_vec=rand(dist)
  # Set up copy of budget
  b=copy(budget.matrix)
  # Loop through species
  for i in eachindex(abun_vec)
    # Get abundance of species
    abun=abun_vec[i]
    # Get species trait
    pref=traits[i]
    # Calculate weighting, giving preference to squares that match with trait
    wv= Vector{Float64}(grid)
    wv[find(reshape(habitat.matrix, (dim[1]*dim[2],1))[grid].==pref)]= 0.9
    wv[find(reshape(habitat.matrix, (dim[1]*dim[2],1))[grid].!=pref)]= 0.1
    # Loop through individuals
      while abun>0
        zs=findin(b[grid], 0)
        deleteat!(grid, zs)
        deleteat!(wv, zs)
        # Randomly choose position on grid (weighted)
      pos=sample(grid, weights(wv))
      # Add individual to this location
      P[i,pos]=P[i,pos]+1
      abun=abun-1
      b[pos]=b[pos]-1
    end
    #sum_abun=mapslices(sum, P, 1)[1,:,:]
    #@rput(sum_abun)
    #R"par(mfrow=c(1,2));image.plot(sum_abun,col=rainbow(50)[1:20], breaks=seq(0,20,1));image(hab, legend = F)"
  end
  # Create MatrixLandscape from P matrix and habitat
  MatrixLandscape(P, habitat, budget)
end


# Function to get the neighbours of a grid square in a matrix in 4 or 8 directions
function get_neighbours(mat::Matrix, y::Int64, x::Int64, chess::Int64=4)
  # Calculate dimensions
  dims=size(mat)
  # Include 4 directions
  if chess==4
    neighbour_vec=[x y-1; x y+1; x-1 y; x+1 y]
  # Include 8 directions
  elseif chess==8
    neighbour_vec=[x y-1; x y+1; x-1 y; x+1 y; x-1 y-1; x-1 y+1; x+1 y-1; x+1 y+1]
  else
    # Give error if other number chosen than 4 or 8
    error("Can only calculate neighbours in 4 or 8 directions")
  end
  # Remove answers outside of the dimensions of the matrix
  remove=vcat(mapslices(all, [neighbour_vec.>=1 neighbour_vec[:,1].<=
    dims[1] neighbour_vec[:,2].<=dims[2]], 2)...)
  neighbour_vec=neighbour_vec[remove,:]
  neighbour_vec
end

# Function to update a Ecosystem after one timestep- stochastic birth, death and movement
function update!(eco::Ecosystem,  birth::Float64, death::Float64, move::Float64,
   l::Float64, s::Float64, timestep::Real)

   # For now keep l>s
   l > s || error("l must be greater than s")
   l >= 0 && s >= 0 || error("l and s must be greater than zero")

  # Calculate abundance in overall grid (to be implemented later?)
  #abun=map(i->sum(eco.partition.abundances[i,:,:]), 1:size(eco.partition.abundances,1))

  # Calculate dimenions of habitat and number of species
  dims = size(eco.partition.abundances)[2:3]
  spp = size(eco.partition.abundances,1)

  # Loop through grid squares
  for x in 1:dims[1]
    for y in 1:dims[2]

      # Get the overall energy budget of that square
      K = eco.abenv.budget.matrix[x, y]
      randomise=collect(1:spp)
      randomise=randomise[randperm(length(randomise))]
      # Loop through species in chosen square
      for j in randomise

        # Get abundances of square we are interested in
        square = eco.partition.abundances[:,x, y]

        if square[j] <= 0
          eco.partition.abundances[j,x, y] = 0
        else
        # Get energy budgets of species in square
        ϵ̄ = eco.spplist.energy.energy
        E = sum(square .* ϵ̄)

        # Alter rates by energy available in current pop & own requirements
        birth_energy = (ϵ̄[j])^(-l-s) * K / E
        death_energy = (ϵ̄[j])^(-l+s) * E / K
        move_energy = 1

        # Calculate effective rates
        birthprob = birth * timestep * birth_energy
        deathprob = death * timestep * death_energy
        moveprob = move * timestep * move_energy

        # If traits are same as habitat type then give birth "boost"
        #if eco.spplist.traits.traits[j] != eco.abenv.habitat.matrix[x, y]
        #  birthrate = birthrate * 0.8
        #end

        # If zero abundance then go extinct
        if square[j] == 0
          birthprob = 0
          deathprob = 0
          moveprob = 0
        end

        # Throw error if rates exceed 1
        #birthprob <= 1 && deathprob <= 1 && moveprob <= 1 ||
        #  error("rates larger than one in binomial draw")


        # Put probabilities into 0 - 1
        probs = map(prob -> 1-exp(-prob), [birthprob, deathprob, moveprob])

        # Calculate how many births and deaths
        births = jbinom(1, Int(square[j]), probs[1])[1]
        deaths = jbinom(1, Int(square[j]), probs[2])[1]

        # Update population
        eco.partition.abundances[j,x, y] = eco.partition.abundances[j,x, y] +
          births - deaths

        # Then calculate movements
        square[j] = eco.partition.abundances[j,x, y]
        moves = jbinom(1, Int(square[j]), probs[3])[1]

        # Update population
        eco.partition.abundances[j,x, y] = eco.partition.abundances[j,x, y] - moves

        # Find neighbours of grid square
        neighbours = get_neighbours(eco.abenv.habitat.matrix, y, x, 8)
        # Randomly sample one of the neighbours
        choose = sample(1:size(neighbours, 1), moves)
        tab = countmap(choose)
        for k in unique(choose)
            destination=neighbours[k, :]
            # Add one to this neighbour
            eco.partition.abundances[j, destination[1], destination[2]] =
              eco.partition.abundances[j, destination[1], destination[2]] + tab[k]
        end
      end
      end
    end
  end
end



# Alternative populate function
function populate!(ml::AbstractStructuredPartition, spplist::SpeciesList,
                   abenv::AbstractAbiotic)
  # Calculate size of habitat
  dim=size(abenv.habitat.matrix)
  grid=collect(1:dim[1]*dim[2])
  # Set up copy of budget
  b=copy(abenv.budget.matrix)
  # Loop through species
  for i in eachindex(spplist.abun)
    # Get abundance of species
    abun=spplist.abun[i]
    # Get species trait
    pref=spplist.traits.traits[i]
    # Calculate weighting, giving preference to squares that match with trait
    wv= Vector{Float64}(grid)
    wv[find(reshape(abenv.habitat.matrix, (dim[1]*dim[2],1))[grid].==pref)]= 0.5
    wv[find(reshape(abenv.habitat.matrix, (dim[1]*dim[2],1))[grid].!=pref)]= 0.5
    # Loop through individuals
      while abun>0
        zs=findin(b[grid], 0)
        deleteat!(grid, zs)
        deleteat!(wv, zs)
        # Randomly choose position on grid (weighted)
      pos=sample(grid, weights(wv))
      # Add individual to this location
      ml.abundances[i,pos]=ml.abundances[i,pos]+1
      abun=abun-1
      b[pos]=b[pos]-1
    end
  end
end

function run_sim(eco::Ecosystem, params::AbstractVector, times::Int64, reps::Int64)

  birth = param[1]
  death = param[2]
  move = param[3]
  timestep = param[4]
  l = param[5]
  s = param[6]

  numSpecies = size(eco.spplist.abun, 1)
  gridSize = size(eco.abenv.habitat.matrix, 1) *  size(eco.abenv.habitat.matrix, 2)
  abun = zeros(times+1, numSpecies, gridSize, reps); ener = zeros(times+1, gridSize, reps)
  for j in 1:reps
    eco_rep = copy_eco(eco)

    for k in 1:gridSize
      abun[1,:, k,  j] = eco_rep.partition.abundances[:, 1]
      ener[1, k,  j] = sum(eco_rep.spplist.abun .* eco_rep.spplist.energy.energy)
    end

    for i in 1:times
        update!(eco_rep, birth, death, move, l, s, timestep)
        for g in 1:gridSize
          abun[i+1, :,g, j] = eco_rep.partition.abundances[: , g]
          ener[i+1, g, j] = sum(eco_rep.partition.abundances[: , g] .* eco_rep.spplist.energy.energy)
        end
    end
    map
  end
  mean_abun = mapslices(mean, abun, 4)
  sd_abun = mapslices(std, abun, 4)
  mean_ener = mapslices(mean, ener, 3)
  sd_ener = mapslices(std, ener, 3)
  [abun, mean_abun, sd_abun, mean_ener, sd_ener]
end

function run_sim(numSpecies::Int64, numTraits::Int64,dist::Distribution,energy_vec,
  numNiches::Int64, grid, K,
   params::AbstractVector,
   times::Int64, reps::Int64)

  birth = param[1]
  death = param[2]
  move = param[3]
  timestep = param[4]
  l = param[5]
  s = param[6]

  gridSize = grid[1] *  grid[2]
  abun = zeros(times+1, numSpecies, gridSize, reps); ener = zeros(times+1, gridSize, reps)

  for j in 1:reps
    sppl = SpeciesList(numSpecies, numTraits, dist,
                       energy_vec)
    abenv = MatrixAbioticEnv(numNiches, grid, K)
    eco = Ecosystem(sppl, abenv)

    for k in 1:gridSize
      abun[1,:, k,  j] = eco.partition.abundances[:, 1]
      ener[1, k,  j] = sum(eco.spplist.abun .* eco.spplist.energy.energy)
    end

    for i in 1:times
        update!(eco, birth, death, move, l, s, timestep)
        for g in 1:gridSize
          abun[i+1, :,g, j] = eco.partition.abundances[: , g]
          ener[i+1, g, j] = sum(eco.partition.abundances[: , g] .* eco.spplist.energy.energy)
        end
    end
    map
  end
  mean_abun = mapslices(mean, abun, 4)
  sd_abun = mapslices(std, abun, 4)
  mean_ener = mapslices(mean, ener, 3)
  sd_ener = mapslices(std, ener, 3)
  [abun, mean_abun, sd_abun, mean_ener, sd_ener]
end
