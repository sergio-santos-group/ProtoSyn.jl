"""
    get_coulomb_potential(d::T; v::Opt{Tuple{T, T, T}} = nothing, kwargs...) where {T <: AbstractFloat}

Return a simple coulomb potential with charges included, as described in
https://www.softschools.com/formulas/physics/potential_energy_electrostatic_point_particles_formula/37/.
The coulomb potential is measured between two charged point-like particles at
distance `d` and charges `kwargs[:qi]` and `kwargs[:qj]`, respectivelly. If a
vector `v` between the two particles is provided, also calculates the resulting
forces.

# See also
[`get_bump_potential`](@ref)

# Examples
```
julia> coulomb = ProtoSyn.Calculators.get_coulomb_potential()

julia> coulomb(2.0, qi = 1.0, qj = -1.0)
-0.5
```
"""
function get_coulomb_potential()
    
    return function coulomb_potential(d::T; v::Opt{Tuple{T, T, T}} = nothing, kwargs...) where {T <: AbstractFloat}
        qi = kwargs[:qi]
        qj = kwargs[:qj]

        e = (qi*qj)/d

        v !== nothing && begin
            _f = (qi*qj)/(d * d)
            f1 = - _f .* v
            f2 =   _f .* v
        end

        v !== nothing && return e, f1, f2
        return e
    end
end