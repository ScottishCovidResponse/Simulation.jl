using RCall
function plot_move(eco::Ecosystem, x::Int64, y::Int64, spp::Int64)
  table = eco.lookup[spp] .+ [x y 0]
  maxGrid = maximum(size(eco.abenv.habitat.matrix))
  # Can't go over maximum dimension
  lower  = find(mapslices(x->all(x.>0), table, 2))
  upper = find(mapslices(x->all(x.<= maxGrid), table, 2))
  valid = intersect(lower, upper)
  table = table[valid, :]
  table[:, 3] = table[:, 3]/sum(table[:, 3])
  table
  A = zeros(size(eco.abenv.habitat.matrix))
  for i in eachindex(table[:, 1])
    A[table[i, 1], table[i, 2]] = table[i, 3]
  end
  @rput A
  R"par(mfrow=c(1,1));library(fields);
  A[A==0]=NA
  image.plot(A)"
end

function plot_abun(abun::AbstractArray, numSpecies::Int64, gridSize::Int64)
  # Plot
  @rput abun
  @rput numSpecies
  @rput gridSize
  R" par(mfrow=c(gridSize,gridSize), mar=c(2, 2, 2, 2))
  for (i in 1:gridSize^2){
      for (k in 1:numSpecies){
        if (k==1) plot_fun=plot else plot_fun=lines
          plot_fun(0:100, abun[, k, 1, i], col=k, xlab='Abundance', ylab='Time', type='l',
          ylim=c(0, max(abun)))
        }
    }"
  end
function plot_divergence(expected::Vector{Float64}, actual::Vector{Float64})
  @rput expected
  @rput actual
  KL = kldivergence(actual, expected); @rput KL
  R"par(mfrow=c(1,1))
    plot(1:length(expected),expected, type='l',
          main = paste('Divergence =', round(KL, 2)), xlab='Abundance',
          ylab='Frequency', ylim=c(0, max(c(expected, actual))))
  abline(h=max(expected), col=1, cex=0.5, lty=3)
  lines(actual, col=2)
  abline(h=max(actual), col=2, cex=0.5, lty=3)
  legend('topright', legend=c('Expected', 'Observed'), col=1:2, pch='-')"
  info("Divergence = ",KL)
end

function plot_divergence(combined::Array{Array{Float64, 1}, 1})
  expected = combined[1]
  actual = combined[2]
  plot_divergence(expected, actual)
end