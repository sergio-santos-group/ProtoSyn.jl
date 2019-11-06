@doc raw"""
    rotation_matrix_from_axis_angle(axis::Vector{Float64}, angle::Float64)

Return a rotation matrix based on the provided axis and angle (in radians).

# Examples
```julia-repl
julia> Aux.rotation_matrix_from_axis_angle([1.1, 2.2, 3.3], π/2)
3×3 Array{Float64,2}:
  0.0714286  -0.658927  0.748808
  0.944641    0.285714  0.16131 
 -0.320237    0.695833  0.642857
```
See also: [`rotate_dihedral!`](@ref Common)
"""
function rotation_matrix_from_axis_angle(axis::Vector{Float64}, angle::Float64)
    # q0 = cos(0.5 * angle)
    # q1, q2, q3 = sin(0.5 * angle) * axis ./ norm(axis)
    λ = sin(0.5 * angle) / norm(axis)
    q0 = cos(0.5 * angle)
    q1 = λ * axis[1]
    q2 = λ * axis[2]
    q3 = λ * axis[3]
    [1-2*q2*q2-2*q3*q3   2*q1*q2-2*q0*q3   2*q1*q3+2*q0*q2;
       2*q2*q1+2*q0*q3 1-2*q3*q3-2*q1*q1   2*q2*q3-2*q0*q1;
       2*q3*q1-2*q0*q2   2*q3*q2+2*q0*q1 1-2*q1*q1-2*q2*q2]
end


function rotation_matrix_from_axis_angle(axis::Vector{Float64}, angle::Float64, R::Array{Float64, 2})
    # q0 = cos(0.5 * angle)
    # q1, q2, q3 = sin(0.5 * angle) * axis ./ norm(axis)
    λ = sin(0.5 * angle) / norm(axis)
    q0 = cos(0.5 * angle)
    q1 = λ * axis[1]
    q2 = λ * axis[2]
    q3 = λ * axis[3]
    R[1,1] = 1-2*q2*q2-2*q3*q3
    R[1,2] =   2*q1*q2-2*q0*q3
    R[1,3] =   2*q1*q3+2*q0*q2
    R[2,1] =   2*q2*q1+2*q0*q3
    R[2,2] = 1-2*q3*q3-2*q1*q1
    R[2,3] =   2*q2*q3-2*q0*q1
    R[3,1] =   2*q3*q1-2*q0*q2
    R[3,2] =   2*q3*q2+2*q0*q1
    R[3,3] = 1-2*q1*q1-2*q2*q2
    R
end


# -------------------------------------------------------------------------------------------------------------------------------------------------------------

@doc raw"""
    calc_dih_angle(a1::Vector{Float64}, a2::Vector{Float64}, a3::Vector{Float64}, a4::Vector{Float64})

Calculates the dihedral angle produced between a1, a2, a3 and a4, in radians.

# Examples
```julia-repl
julia> Aux.calc_dih_angle([1.0, 1.0, 1.0], [2.1, 2.1, 2.1], [3.0, 2.0, 5.0], [5.0, 5.0, 5.0])
3.141592653589793
```
See also: [`apply_initial_conf!`](@ref Common)
"""
function calc_dih_angle(a1::Vector{Float64}, a2::Vector{Float64}, a3::Vector{Float64}, a4::Vector{Float64})::Float64

    b1 = a2 - a1
    b2 = a3 - a2
    b3 = a4 - a3
    n1 = cross(b1, b2)
    n2 = cross(b2, b3)
    x = dot(cross(n1, n2), b2) / sqrt(dot(b2, b2))
    y = dot(n1, n2)
    return atan(x, y)
end

function calc_angle(a1::Vector{Float64}, a2::Vector{Float64}, a3::Vector{Float64})::Float64

    v21 = a1 - a2
    v23 = a3 - a2
    return acos(dot(v21, v23) / (norm(v21) * norm(v23)))
end


@doc raw"""
    rand_vector_in_sphere()::Vector{Float64}

Return a random vector in a sphere.

# Examples
```julia-repl
julia> Aux.rand_vector_in_sphere()
3-element Array{Float64,1}:
  0.34470440733141694
 -0.746086857969672  
  0.5696782178837796
```
"""
function rand_vector_in_sphere()::Vector{Float64}
    theta::Float64 = 2 * π * rand()
    phi::Float64 = acos(1 - 2 * rand())
    x::Float64 = sin(phi) * cos(theta)
    y::Float64 = sin(phi) * sin(theta)
    z::Float64 = cos(phi)
    return [x, y, z]
end

@doc raw"""
    linreg(x::Union{Vector{Float64}, Vector{Int64}}, y::Union{Vector{Float64}, Vector{Int64}})::Float64

Returns the linear regression slope value for the x and y values.

# Examples
```julia-repl
julia> Aux.linreg([1, 2, 3], [1, 2, 3])
1.0
```
"""
function linreg(x::Union{Vector{Float64}, Vector{Int64}}, y::Union{Vector{Float64}, Vector{Int64}})::Float64

    x = convert(Vector{Float64}, x)
    y = convert(Vector{Float64}, y)

    avgx = sum(x)/length(x)
    avgy = sum(y)/length(y)

    x_avgx = x .- avgx
    y_avgy = y .- avgy
    x_avgx_y_avgy = x_avgx .* y_avgy
    x_avgx_sq = x_avgx .* x_avgx
    y_avgy_sq = y_avgy .* y_avgy

    r = sum(x_avgx_y_avgy) / sqrt(sum(x_avgx_sq) * sum(y_avgy_sq))
    sy = sqrt(sum(y_avgy_sq) / (length(y) - 1))
    sx = sqrt(sum(x_avgx_sq) / (length(x) - 1))

    return r * (sy/sx)
end


@inline function get_max_force(n_atoms::Int64, mat::Array{T, 2}) where {T<:AbstractFloat}
    max_force = zero(eltype(mat))
    for i=1:n_atoms
        forceSq = mat[i, 1] ^ 2 + mat[i, 2] ^ 2 + mat[i, 3] ^ 2
        if forceSq > max_force
            max_force = forceSq
        end
    end
    return sqrt(max_force)
end