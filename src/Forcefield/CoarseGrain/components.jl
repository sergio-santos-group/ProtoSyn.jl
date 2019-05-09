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
mutable struct SolvPair <: Abstract.ForcefieldComponent
    i::Int64
    coef::Float64
end
Base.show(io::IO, b::SolvPair) = print(io, "Forcefield.CoarseGrain.SolvPair(i=$(b.i), coef=$(b.coef))")


@doc raw"""
    HbPair(donors::Vector{HbPair}, acceptors::Vector{HbPair}, coef::Float64)

Hydrogen Bonding pair.

# Arguments
- `charged::Int64`: *global* index for the charged atom in the HbPair (ex: hydrogen in N-H or oxygen in C=O pairs).
- `base::Int64`: *global* index for the anchor atom in the HbPair (ex: nitrogen in N-H or carbon in C=O pairs).

# Examples
```julia-repl
julia> Forcefield.CoarseGrain.HbPair(13, 14)
Forcefield.CoarseGrain.HbGroup(charged=13, base=14)
```
See algo: [`CoarseGrain.evaluate!`](@ref)
"""
mutable struct HbPair
    charged::Int64
    base::Int64
end
Base.show(io::IO, b::HbPair) = print(io, "Forcefield.CoarseGrain.HbPair(charged=$(b.charged), base=$(b.base))")


@doc raw"""
    HbNetwork(donors::Vector{HbPair}, acceptors::Vector{HbPair}, coef::Float64)

Hydrogen Bonding network, defined by the combination os all `donors` and `acceptors` here defined. The energy contribution of this
`HbNetwork` is multiplied by the defined `coef`.

# Arguments
- `donors::Vector{HbPair}`: List of hydrogen bond donor pairs.
- `acceptors::Vector{HbPair}`: List of hydrogen bond acceptor pairs.
- `coef::Float64`: λ hydrogen bond energy coeficient.

# Examples
```julia-repl
julia> Forcefield.CoarseGrain.HbNetwork(donors, acceptors, 1.0)
Forcefield.CoarseGrain.HbGroup(donors=10, acceptors=10, coef=1.0)
```
See algo: [`CoarseGrain.evaluate!`](@ref)
"""
mutable struct HbNetwork <: Abstract.ForcefieldComponent
    donors::Vector{HbPair}
    acceptors::Vector{HbPair}
    coef::Float64
end
HbNetwork(; coef = 0.0) = HbNetwork(Vector{HbPair}(), Vector{HbPair}(), coef)
Base.show(io::IO, b::HbNetwork) = print(io, "Forcefield.CoarseGrain.HbNetwork(donors=$(length(b.donors)), acceptors=$(length(b.acceptors)), coef=$(b.coef))")
Base.length(b::HbNetwork) = (length(b.donors), length(b.acceptors))