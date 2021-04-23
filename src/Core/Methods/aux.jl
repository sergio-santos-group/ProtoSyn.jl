using Statistics

"""
    rand_vector_in_sphere([::Type{T}]) where {T <: AbstractFloat}

Return a random vector in a sphere, with norm = 1. The X, Y and Z cartesian
coordinates are of type `T`. If not provided, the `ProtoSyn.Units.defaultFloat`
type will be used.

# Examples
```jldoctest
julia> ProtoSyn.rand_vector_in_sphere()
3-element Array{Float64,1}:
 0.4182585394795512
 0.4965570141814753
 0.7605859095589471
```
"""
function rand_vector_in_sphere(::Type{T}) where {T <: AbstractFloat}
    theta = 2 * π * rand(T)
    phi   = acos(1 - 2 * rand(T))
    x     = sin(phi) * cos(theta)
    y     = sin(phi) * sin(theta)
    z     = cos(phi)
    return Vector{T}([x, y, z])
end

rand_vector_in_sphere() = rand_vector_in_sphere(ProtoSyn.Units.defaultFloat)


"""
    rotation_matrix_from_axis_angle(axis::Vector{T}, angle::T) where {T <: AbstractFloat}

Return a rotation matrix based on the provided axis and angle (in radians).

# Examples
```jldoctest
julia rotation_matrix_from_axis_angle([1.1, 2.2, 3.3], π/2)
3×3 Array{Float64, 2}:
  0.0714286  -0.658927  0.748808
  0.944641    0.285714  0.16131 
 -0.320237    0.695833  0.642857
```
"""
function rotation_matrix_from_axis_angle(axis::Vector{T}, angle::T) where {T <: AbstractFloat}
    q0 = cos(0.5 * angle)
    q1, q2, q3 = sin(0.5 * angle) * axis ./ norm(axis)
    return [1-2*q2*q2-2*q3*q3   2*q1*q2-2*q0*q3   2*q1*q3+2*q0*q2;
            2*q2*q1+2*q0*q3   1-2*q3*q3-2*q1*q1   2*q2*q3-2*q0*q1;
            2*q3*q1-2*q0*q2     2*q3*q2+2*q0*q1 1-2*q1*q1-2*q2*q2]
end


"""
    unit_circle(value::T) where {T <: AbstractFloat}

Maps an angle value (in radians) from ]-∞, +∞[ to [0, 2π[ range. Returns value
in radians.

# Examples
```jldoctest
julia ProtoSyn.unit_circle(-2pi)
0.0

julia> ProtoSyn.unit_circle(-(3/2)pi)
1.5707963267948966
```
"""
function unit_circle(value::T) where {T <: AbstractFloat}
    r = deg2rad((value > 0 ? value : (2*pi + value)) * 360 / (2*pi))
    return mod(r, 2*pi)
end


"""
    gpu_allocation()

Return the current fraction of the GPU memory allocated (in range [0, 1]).

# Examples
```jldoctest
julia> ProtoSyn.gpu_allocation()
0.07821749405343822
```
"""
function gpu_allocation()
    mem = (CUDA.total_memory() - CUDA.available_memory()) / CUDA.total_memory()
    return ProtoSyn.Units.defaultFloat(mem)
end