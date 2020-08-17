using Distributions
"""
    ExponentialDecay{TR} <: AbstractTraitRelationship{TR}

The relationship between a virus trait and its environment, paramaterised on any TR.
"""
mutable struct ExponentialDecay{TR} <: AbstractTraitRelationship{TR}
end

function (::ExponentialDecay{TR})(current::TR, opt::TR) where TR
    return exp((current - opt)/unit(current))
end
iscontinuous(tr::ExponentialDecay{TR}) where TR = true
function eltype(tr::ExponentialDecay{TR}) where TR
    return TR
end

"""
    ExponentialTrait{C <: Number} <: ContinuousTrait{C}

Trait type that holds Gaussian mean and variance trait information for each species, of any number type `C`.
"""
mutable struct ExponentialTrait{C <: Number} <: ContinuousTrait{C}
  opt::Array{C, 1}
end
iscontinuous(trait::ExponentialTrait{C}) where C = true
function eltype(trait::ExponentialTrait{C}) where C
    return C
end
function ExponentialTrait(opt::Array{C, 1}) where C  <: Unitful.Temperature
    optK = uconvert.(K, opt)
    return ExponentialTrait{typeof(1.0K)}(optK)
end
