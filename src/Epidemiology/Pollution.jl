
abstract type AbstractPollution end

mutable struct NoPollution <: AbstractPollution end

function _get_pollution(pol::NoPollution, j::Int64)
    return 1.0 * μg * m^-3
end

mutable struct GriddedPollution <: AbstractPollution
    matrix::AxisArray{typeof(1.0μg*m^-3), 2}
end

function _get_pollution(pol::GriddedPollution, j::Int64)
    width = size(pol.matrix, 1)
    x, y = convert_coords(j, width)
    return pol.matrix[x, y]
end
