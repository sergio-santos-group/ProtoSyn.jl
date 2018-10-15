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
See also: [`Mutators.Dihedral.rotate_dihedral!`](@ref Mutators)
"""
function rotation_matrix_from_axis_angle(axis::Vector{Float64}, angle::Float64)
    q0 = cos(0.5 * angle)
    q1, q2, q3 = sin(0.5 * angle) * axis ./ norm(axis)
    [1-2*q2*q2-2*q3*q3   2*q1*q2-2*q0*q3   2*q1*q3+2*q0*q2;
       2*q2*q1+2*q0*q3 1-2*q3*q3-2*q1*q1   2*q2*q3-2*q0*q1;
       2*q3*q1-2*q0*q2   2*q3*q2+2*q0*q1 1-2*q1*q1-2*q2*q2]
end

# TODO: Document function
function calc_dih_angle(a1::Vector{Float64}, a2::Vector{Float64}, a3::Vector{Float64}, a4::Vector{Float64};
    radians::Bool = false)

    b1 = a2 - a1
    b2 = a3 - a2
    b3 = a4 - a3
    n1 = cross(b1, b2)
    n2 = cross(b2, b3)
    x = dot(cross(n1, n2), b2) / sqrt(dot(b2, b2))
    y = dot(n1, n2)
    v = atan(x, y)

    if radians
        return v
    else
        return rad2deg(v)
    end
end