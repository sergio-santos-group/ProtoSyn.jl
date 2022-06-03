"""
    get_harmonic_potential(a::T, b::T, c::T) where {T <: AbstractFloat}

Return a simple harmonic potential with the following equation:

> U = ((d+a)²-b)*c

The harmonic potential is measured between two point-like particles at distance
`d` (doesn't take charges into consideration). If a vector `v` between the two
particles is provided, also calculates the resulting forces.

# See also
[`get_bump_potential`](@ref)

# Examples
```
julia> harmonic = ProtoSyn.Calculators.get_harmonic_potential(1.0, 1.0, 1.0)

julia> harmonic(1.0)
3.0
```
"""
function get_harmonic_potential(a::T, b::T, c::T) where {T <: AbstractFloat}
    
    return function harmonic_potential(d::T; v::Opt{Tuple{T, T, T}} = nothing, kwargs...) where {T <: AbstractFloat}

        e = (((d+a)^2)-b)*c

        v !== nothing && begin
            _f = 2*c*(x-a)
            f1 = - _f .* v
            f2 =   _f .* v
        end

        v !== nothing && return e, f1, f2
        return e
    end
end