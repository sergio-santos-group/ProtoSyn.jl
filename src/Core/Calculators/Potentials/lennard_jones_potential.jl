"""
    get_lennard_jones_potential(d::T; v::Opt{Tuple{T, T, T}} = nothing, kwargs...) where {T <: AbstractFloat}

Return a simple ennard-Jones potential. The coulomb potential is measured
between two non-charged point-like particles at distance `d`. If a vector `v`
between the two particles is provided, also calculates the resulting forces.

# See also
[`get_bump_potential`](@ref)

# Examples
```
julia> lj = ProtoSyn.Calculators.get_lennard_jones_potential()

julia> lj(2.0)
-0.015380859375
```
"""
function get_lennard_jones_potential()
    
    return function lennard_jones_potential(d::T; v::Opt{Tuple{T, T, T}} = nothing, kwargs...) where {T <: AbstractFloat}

        dij_sq = d * d
        lj2    = T(1) / dij_sq
        lj6    = lj2*lj2*lj2
        e      = T(1)*(lj6*lj6 - lj6)

        v !== nothing && begin
            _f = T(24) * (lj6 - T(2) * lj6 * lj6) / dij_sq
            f1 = - _f .* v
            f2 =   _f .* v
        end

        v !== nothing && return e, f1, f2
        return e
    end
end