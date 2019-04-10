module Sidechain

using StatsBase
using Random
using ..Common
using ..Aux

rng = MersenneTwister(1234)

@doc raw"""
    SidechainMutator(sidechains::Vector{Common.SidechainMetadata}, p_mut::Float64)

Holds all necessary parameters for the correct simulation of sidechain movements.

# Arguments
- `sidechains::Vector{Common.SidechainMetadata}`: List of sidechains avaliable to be mutated in sidechain movements.
- `p_mut::Float64`: Probability of mutation of each sidechain.

# Examples
```julia-repl
julia> Mutators.Sidechain.SidechainMutator(sidechains, 0.25)
SidechainMutator(sidechains=68, p_pmut=0.05)
```
See also: [`run!`](@ref)
"""
mutable struct SidechainMutator
    sidechains::Vector{Common.SidechainMetadata}
    p_mut::Float64
end
SidechainMutator(sidechains::Vector{Common.SidechainMetadata}, p_mut=0.0) = SidechainMutator(sidechains, p_mut)
Base.show(io::IO, b::SidechainMutator) = print(io, "SidechainMutator(sidechains=$(length(b.sidechains)), p_mut=$(b.p_mut)")


@doc raw"""
    run!(state::Common.State, mutator::SidechainMutator)

Iterate over a list of [`Common.SidechainMetadata`](@ref) and perform sidechain mutations on the current
[`Common.State`](@ref). The probability of each sidechain undergo movements is defined in the
[`SidechainMutator`](@ref).`p_mut`. The new list of angles (one of each `chi` in the sidechain) is obtained
from [`SidechainMetadata`](@ref).`rotamers`, and sampled with a random weighted probability, defined in
[`SidechainMetadata`](@ref).`weights`. After mutation, the [`Common.State`](@ref) is updated with the new conformation.
Returns the number of mutations performed.

# Examples
```julia-repl
julia> Mutators.Sidechain.run!(state, mutator)
```
See also: [`Common.rotate_dihedral_to!`](@ref Common)
"""
@inline function run!(state::Common.State, mutator::SidechainMutator)
    
    count::Int64 = 0
    for sidechain in mutator.sidechains
        if rand() < mutator.p_mut
            rotamer = sidechain.rotamers[sample(1:length(sidechain.weights), Weights(sidechain.weights))]
            for chi_index in 1:length(rotamer.chis)
                chi_center = rotamer.chis[chi_index]
                chi_range  = rotamer.ranges[chi_index]
                chi_value  = (chi_range/1.175) * randn(rng, Float64, 1) .+ chi_center
                Common.rotate_dihedral_to!(state.xyz, sidechain.dihedrals[chi_index], chi_value[1])
            end
            count += 1
        end
    end
    return count
end

end