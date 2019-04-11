@doc raw"""
    SolvPair(i::Int64, coef::Float64)

Solvation pair.

# Arguments
- `i::Int64`: *global* atom index of the residue Cα;
- `coef::Float64`: aminoacid solvation coeficient.

# Examples
```julia-repl
julia> Forcefield.CoarseGrain.SolvPair(1, -3.5)
Forcefield.CoarseGrain.SolvPair(i=1, coef=-3.5)
```
See algo: [`CoarseGrain.evaluate!`](@ref)
"""
mutable struct SolvPair
    i::Int64
    coef::Float64
end
Base.show(io::IO, b::SolvPair) = print(io, "Forcefield.CoarseGrain.SolvPair(i=$(b.i), coef=$(b.coef))")


@doc raw"""
    HbGroup(n::Int64, h::Int64, c::Int64, o::Int64)

Hydrogen Bonding group.

# Arguments
- `n::Int64, h::Int64, c::Int64, o::Int64`: *global* indices of the N, H, C and O atoms of one residue;
- `coef::Float64`: aminoacid hydrogen bonding coeficient.

# Examples
```julia-repl
julia> Forcefield.CoarseGrain.HbGroup(1, 2, 4, 5, 1.0)
Forcefield.CoarseGrain.HbGroup(n=1, h=2, c=4, o=5, coef=1.0)
```
See algo: [`CoarseGrain.evaluate!`](@ref)
"""
mutable struct HbGroup
    n::Int64
    h::Int64
    c::Int64
    o::Int64
    coef::Float64
end
Base.show(io::IO, b::HbGroup) = print(io, "Forcefield.CoarseGrain.HbGroup(N=$(b.n), H=$(b.h), C=$(b.c), O=$(b.o), coef=$(b.coef))")
