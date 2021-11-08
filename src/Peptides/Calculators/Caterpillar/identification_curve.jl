"""
# TODO
"""
function linear(distance::T; rmax::T = 9.0, sc::T = 1.0) where {T <: AbstractFloat}
    return distance > rmax ? 0.0 : - (distance/rmax) + 1
end


"""
# TODO
"""
function sigmoid(distance::T; rmax::T = 9.0, sc::T = 1.0) where {T <: AbstractFloat}
    return 1 - (1 / (1 + exp(sc * (rmax - distance))))
end


"""
# TODO
"""
function sigmoid_normalized(distance::T; rmax::T = 9.0, sc::T = 1.0) where {T <: AbstractFloat}
    return (1 - (1 / (1 + exp(sc * (rmax - distance)))))/distance
end