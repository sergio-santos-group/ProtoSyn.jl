@doc raw"""
    SolvPair(i::Int64, coef::Float64)

Solvation pair.

# Arguments
- `i::Int64`: *global* atom index of the residue Cα.
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