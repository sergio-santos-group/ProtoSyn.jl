"""
    get_flat_bottom_potential(;d1::T = 0.0, d2::T = 0.0, d3::T = Inf, d4::T = Inf) where {T <: AbstractFloat}

Return a flat-bottom potential function, using the specified distances. The
potential is made up of 5 different sectors, each with the following functions:

\$f_{1}) \\;\\;\\;\\;\\;\\; e = m_{1} \\cdot d + b_{1} \\;\\;\\;\\;\\;\\; \\left \\{ d < d_{1} \\right \\}\\,\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\$

\$f_{2}) \\;\\;\\;\\;\\;\\; e = \\left (d-d_{2}  \\right )^{2} \\;\\;\\;\\;\\;\\;\\;\\; \\left \\{ d_{1} \\leqslant d < d_{2} \\right \\}\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\$

\$f_{3}) \\;\\;\\;\\;\\;\\; e = 0 \\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\; \\left \\{ d_{2} \\leqslant d \\leqslant d_{3} \\right \\}\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\$

\$f_{4}) \\;\\;\\;\\;\\;\\; e = \\left (d-d_{3}  \\right )^{2} \\;\\;\\;\\;\\;\\;\\;\\; \\left \\{ d_{4} < d \\leqslant d_{4} \\right \\}\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\$

\$f_{5}) \\;\\;\\;\\;\\;\\; e = m_{2} \\cdot d + b_{2} \\;\\;\\;\\;\\;\\;  \\left \\{ d > d_{4} \\right \\}\\,\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\$

Where 

\$m_{1} = 2 \\left ( d_{1}-d_{2} \\right ) \\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\$
\$b_{1} = f_{2}\\left ( d_{1} \\right ) - m_{1} \\cdot d_{1} \\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\$
\$m_{2} = 2\\left ( d_{4} - d_{3} \\right ) \\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\$
\$b_{2} = f_{4}\\left ( d_{4} \\right ) - m_{2} \\cdot d_{4} \\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\$

*The resulting function can be called with the following signature:*

```
flat_bottom_potential(d::T; v::Opt{Vector{T}} = nothing) where {T <: AbstractFloat}
```

Return an energy value based on the provided distance `d`. If a vector `v` is
also provided (optional), the flat-bottom restraint will also return the forces
`f1` and `f2` (the forces felt on both ends of the vector `v`). The vector `v`
should have length = 3, corresponding to the 3 dimensions of the distance
between the two [`Atom`](@ref) instances (X, Y and Z). For more information on
the flat-bottom potential, please read:
[https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4692055/](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4692055/). 

# See also
[`apply_potential!`](@ref ProtoSyn.Calculators.apply_potential!) [`calc_flat_bottom_restraint`](@ref ProtoSyn.Calculators.Restraints.calc_flat_bottom_restraint)

# Examples
```jldoctest
julia> f = ProtoSyn.Calculators.get_flat_bottom_potential(d1 = 1.0, d2 = 2.0, d3 = 3.0, d4 = 4.0);

julia> f(2.5)
0.0

julia> f(1.73, v = (1.0, 1.0, 1.0))
(0.0729, (0.54, 0.54, 0.54), (-0.54, -0.54, -0.54))
```
"""
function get_flat_bottom_potential(A::Union{Type{ProtoSyn.SISD_0}, Type{ProtoSyn.CUDA_2}}; d1::T = 0.0, d2::T = 0.0, d3::T = Inf, d4::T = Inf) where {T <: AbstractFloat}

    @assert d2 >= d1 "d2 must be of equal value or over d1"
    @assert d3 >= d2 "d3 must be of equal value or over d2"
    @assert d4 >= d3 "d4 must be of equal value or over d3"
    @assert d1 !== Inf && d2 !== Inf "Neither `d1` or `d2` can be infinite"

    return function flat_bottom_potential(d::T; v::Opt{Tuple{T, T, T}} = nothing, qi::T = 0.0, qj::T = 0.0) where {T <: AbstractFloat}

        # Note: 'qi' and 'qj' refer to the atomic partial charges, which are not
        # necessary for this potential. They are included for uniformization
        # purposes only.
        
        v !== nothing && begin
            f1 = (T(0.0), T(0.0), T(0.0))
            f2 = (T(0.0), T(0.0), T(0.0))
        end

        
        d12     = d1 - d2
        m1      = 2 * d12
        b1      = d12 * d12 - m1 * d1
        
        d43     = d4 - d3
        m2      = 2 * d43
        b2      = d43 * d43 - m2 * d4
        
        if d < d1 # * Linear Left (1)
            e = m1 * d + b1
            
            v !== nothing && begin
                f1 =    v .* m1
                f2 = .- v .* m1
            end
        elseif d1 <= d < d2 # * Quadratic Left (2)
            e = (d - d2) * (d - d2)
            
            v !== nothing && begin
                factor1 = 2 * (d - d2)
                f1      =    v .* factor1
                f2      = .- v .* factor1
            end
        elseif d2 <= d <= d3 # * Flat (3)
            e = T(0)
        elseif d3 < d <= d4 # * Quadratic Right (4)
            e = (d - d3) * (d - d3)
                
            v !== nothing && begin
                factor2 = 2 * (d - d3)
                f1      =    v .* factor2
                f2      = .- v .* factor2
            end            
        else # * Linear Right (5)
            e = m2 * d + b2
            
            v !== nothing && begin
                f1 =    v .* m2
                f2 = .- v .* m2
            end
        end

        v !== nothing && return e, f1, f2
        return e
    end
end


function get_flat_bottom_potential(A::Type{ProtoSyn.SIMD_1}; d1::T = 0.0, d2::T = 0.0, d3::T = Inf, d4::T = Inf) where {T <: AbstractFloat}

    @assert d2 >= d1 "d2 must be of equal value or over d1"
    @assert d3 >= d2 "d3 must be of equal value or over d2"
    @assert d4 >= d3 "d4 must be of equal value or over d3"
    @assert d1 !== Inf && d2 !== Inf "Neither `d1` or `d2` can be infinite"

    return function flat_bottom_potential(d::T; v::Opt{Vec{4, T}} = nothing, qi::T = 0.0, qj::T = 0.0) where {T <: AbstractFloat}

        # Note: 'qi' and 'qj' refer to the atomic partial charges, which are not
        # necessary for this potential. They are included for uniformization
        # purposes only.

        v !== nothing && begin
            f1 = Vec((T(0.0), T(0.0), T(0.0)))
            f2 = Vec((T(0.0), T(0.0), T(0.0)))
        end

        
        d12     = d1 - d2
        m1      = 2 * d12
        b1      = d12 * d12 - m1 * d1
        
        d43     = d4 - d3
        m2      = 2 * d43
        b2      = d43 * d43 - m2 * d4
        
        if d < d1 # * Linear Left (1)
            e = m1 * d + b1
            
            v !== nothing && begin
                f1 =   v * m1
                f2 = - v * m1
            end
        elseif d1 <= d < d2 # * Quadratic Left (2)
            e = (d - d2) * (d - d2)
            
            v !== nothing && begin
                factor1 = 2 * (d - d2)
                f1      =   v * factor1
                f2      = - v * factor1
            end
            elseif d2 <= d <= d3 # * Flat (3)
                e = T(0)
            elseif d3 < d <= d4 # * Quadratic Right (4)
                e = (d - d3) * (d - d3)
                
                v !== nothing && begin
                factor2 = 2 * (d - d3)
                f1      =   v * factor2
                f2      = - v * factor2
            end            
        else # * Linear Right (5)
            e = m2 * d + b2
            
            v !== nothing && begin
            f1 =   v * m2
            f2 = - v * m2
        end
    end

    v !== nothing && return e, f1, f2
        return e
    end
end

get_flat_bottom_potential(; d1::T = 0.0, d2::T = 0.0, d3::T = Inf, d4::T = Inf) where {T <: AbstractFloat} = begin
    get_flat_bottom_potential(ProtoSyn.acceleration.active; d1 = d1, d2 = d2, d3 = d3, d4 = d4)
end