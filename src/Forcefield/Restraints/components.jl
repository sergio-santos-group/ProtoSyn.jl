@doc raw"""
    DistanceFBR(a1::Int64, a2::Int64, r1::Float64, r2::Float64, r3::Float64, r4::Float64, c::Float64)

Distance flat-bottomed restraint

# Arguments
- `a1::Int64, a2::Int64`: *global* atom indices.
- `r1::Float64, r2::Float64, r3::Float64, r4::Float64`: Distances that set the 5 zones of a classical flat-bottomed restraint.
- `c::Float64`: force constant (kJ mol⁻¹ nm⁻²).

# Examples
```julia-repl
julia> Forcefield.HarmonicAngle(1, 2, 3, 670.0, 1.92)
Forcefield.HarmonicAngle(a1=1, a2=2, a3=3, k=670.0, θ=1.92)
```
See algo: [`Amber.evaluate!`](@ref Forcefield)
"""
struct DistanceFBR
    a1::Int64   # in nm
    a2::Int64   # in nm
    r1::Float64 # in rad
    r2::Float64 # in rad
    r3::Float64 # in rad
    r4::Float64 # in rad
    c::Float64
end
Base.show(io::IO, b::DistanceFBR) = print(io, "Forcefield.Restraint.DistanceFBR(a1=$(b.a1), a2=$(b.a2), r1=$(b.r1), r2=$(b.r2), r3=$(b.r3), r4=$(b.r4), c=$(b.c))")


struct DihedralFBR
    a1::Int64   # in nm
    a2::Int64   # in nm
    a3::Int64   # in nm
    a4::Int64   # in nm
    r1::Float64 # in rad
    r2::Float64 # in rad
    r3::Float64 # in rad
    r4::Float64 # in rad
    c::Float64
end
Base.show(io::IO, b::DihedralFBR) = print(io, "Forcefield.Restraint.DihedralFBR(a1=$(b.a1), a2=$(b.a2), a3=$(b.a3), a4=$(b.a4), r1=$(rad2deg(b.r1)), r2=$(rad2deg(b.r2)), r3=$(rad2deg(b.r3)), r4=$(rad2deg(b.r4)), c=$(b.c))")