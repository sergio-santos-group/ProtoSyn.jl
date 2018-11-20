@doc raw"""
    DistanceFBR(a1::Int64, a2::Int64, r1::Float64, r2::Float64, r3::Float64, r4::Float64, c::Float64)

Distance flat-bottomed restraint

# Arguments
- `a1::Int64, a2::Int64`: *global* atom indices.
- `r1::Float64, r2::Float64, r3::Float64, r4::Float64`: Distances that set the 5 zones of a classical flat-bottomed restraint (in nm).
- `c::Float64`: force constant (kJ mol⁻¹ nm⁻²).

# Examples
```julia-repl
julia> Forcefield.Restraints.DistanceFBR(1, 2, 2.0, 4.0, 6.0, 8.0, 1e4)
Forcefield.Restraints.DistanceFBR(a1=1, a2=2, r1=2.0, r2=4.0, r3=6.0, r4=8.0, c=1e4)
```
See algo: [`Restraint.evaluate!`](@ref Forcefield)
"""
struct DistanceFBR
    a1::Int64   # global index
    a2::Int64   # global index
    r1::Float64 # in nm
    r2::Float64 # in nm
    r3::Float64 # in nm
    r4::Float64 # in nm
    c::Float64
end
Base.show(io::IO, b::DistanceFBR) = print(io, "Forcefield.Restraints.DistanceFBR(a1=$(b.a1), a2=$(b.a2), r1=$(b.r1), r2=$(b.r2), r3=$(b.r3), r4=$(b.r4), c=$(b.c))")


@doc raw"""
    DihedralFBR(a1::Int64, a2::Int64, a3::Int64, a4::Int64, r1::Float64, r2::Float64, r3::Float64, r4::Float64, c::Float64)

Distance flat-bottomed restraint

# Arguments
- `a1::Int64, a2::Int64, a3::Int64, a4::Int64`: *global* atom indices.
- `r1::Float64, r2::Float64, r3::Float64, r4::Float64`: Angles that set the 5 zones of a classical flat-bottomed restraint (in rad).
- `c::Float64`: force constant (kJ mol⁻¹ rad⁻²).

# Examples
```julia-repl
julia> Forcefield.Restraints.DihedralFBR(1, 2, 3, 4, -π, -π/2, π/2, π, 1e4)
Forcefield.Restraints.DistanceFBR(a1=1, a2=2, a3=3, a4=4, r1=-180, r2=-90, r3=90, r4=180, c=1e4)
```
See algo: [`Restraint.evaluate!`](@ref Forcefield)
"""
struct DihedralFBR
    a1::Int64   # global index
    a2::Int64   # global index
    a3::Int64   # global index
    a4::Int64   # global index
    r1::Float64 # in rad
    r2::Float64 # in rad
    r3::Float64 # in rad
    r4::Float64 # in rad
    c::Float64

    function DihedralFBR(a1::Int64, a2::Int64, a3::Int64, a4::Int64, r1::Float64, r2::Float64, r3::Float64, r4::Float64, c::Float64)
        if !(r1 < r2 < r3 < r4)
            error("The DihedralFBR angles must obey the following sequence: r1 < r2 < r3 < r4.") 
        end
        new(dihedrals, angle_sampler, p_mut, step_size)
    end
end
Base.show(io::IO, b::DihedralFBR) = print(io, "Forcefield.Restraints.DihedralFBR(a1=$(b.a1), a2=$(b.a2), a3=$(b.a3), a4=$(b.a4), r1=$(rad2deg(b.r1)), r2=$(rad2deg(b.r2)), r3=$(rad2deg(b.r3)), r4=$(rad2deg(b.r4)), c=$(b.c))")