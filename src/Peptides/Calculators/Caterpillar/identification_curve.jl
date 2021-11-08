"""
    linear(distance::T; rmax::T = 9.0, sc::T = 1.0) where {T <: AbstractFloat}

Returns the linear identification curve weight `w1` (max at `distance` 0.0,
linearly decreasing until `rmax`). Note that slope control `sc` has no effect.

# See also
[`sigmoid`](@ref) [`sigmoid_normalized`](@ref)

# Examples
```
julia> ProtoSyn.Peptides.Calculators.Caterpillar.linear(10.0, rmax = 20.0)
0.5
```
"""
function linear(distance::T; rmax::T = 9.0, sc::T = 1.0) where {T <: AbstractFloat}
    return distance > rmax ? 0.0 : - (distance/rmax) + 1
end


"""
    sigmoid(distance::T; rmax::T = 9.0, sc::T = 1.0) where {T <: AbstractFloat}

Returns the sigmoid identification curve weight `w1` (max at `distance` 0.0,
min at distance `rmax`). The slope control `sc` defines how sharp the sigmoid
curve is.

# See also
[`linear`](@ref) [`sigmoid_normalized`](@ref)

# Examples
```
julia> ProtoSyn.Peptides.Calculators.Caterpillar.sigmoid(10.0, rmax = 20.0)
0.9999546021312976
```
"""
function sigmoid(distance::T; rmax::T = 9.0, sc::T = 1.0) where {T <: AbstractFloat}
    return 1 - (1 / (1 + exp(sc * (rmax - distance))))
end


"""
    sigmoid_normalized(distance::T; rmax::T = 9.0, sc::T = 1.0) where {T <: AbstractFloat}

Returns the sigmoid identification curve weight `w1` (max at `distance` 0.0,
min at distance `rmax`), normalized by the given distance (should only be 
applied in Neighbor Vector (NV) algorithms). The slope control `sc` defines how
sharp the sigmoid curve is.

# See also
[`linear`](@ref) [`sigmoid`](@ref)

# Examples
```
julia> ProtoSyn.Peptides.Calculators.Caterpillar.sigmoid_normalized(10.0, rmax = 20.0)
0.09999546021312976
```
"""
function sigmoid_normalized(distance::T; rmax::T = 9.0, sc::T = 1.0) where {T <: AbstractFloat}
    return (1 - (1 / (1 + exp(sc * (rmax - distance)))))/distance
end