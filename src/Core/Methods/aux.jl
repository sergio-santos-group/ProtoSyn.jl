using Statistics

"""
    rand_vector_in_sphere(::Type{T}) where {T <: AbstractFloat}

Return a random vector in a sphere, with norm = 1.

# Examples
```julia-repl
julia> rand_vector_in_sphere()
3-element Array{Float64, 1}:
  0.34470440733141694
 -0.746086857969672  
  0.5696782178837796
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
```julia-repl
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

Maps an angle value (in radians) from ]-∞, +∞[ to [0, 2π] range. Returns value
in radians.

# Examples
```julia-repl
julia unit_circle(-2pi)
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
"""
function center_of_mass(pose::Pose, idxs::Vector{Int})
    return mean(pose.state.x[:, idxs], dims = 2)
end

function random_atom(pose::Pose, idxs::Vector{Int})
    return pose.state.x[:, rand(idxs)]
end