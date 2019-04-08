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

Iterate over a list of [`Common.Dihedral`](@ref) (`dihedrals`) and perform dihedral movements on the current
[`Common.State`](@ref). The probability of each dihedral undergo movements is defined in the
[`DihedralMutator`](@ref).`p_mut`. The new angle is obtained from [`DihedralMutator`](@ref).`angle_sampler`, who should
return a `Float64` in radians.
After movement, the [`Common.State`](@ref) is updated with the new conformation.
Returns the number of rotations performed.

# Examples
```julia-repl
julia> Mutators.Dihedral.run!(state, mutator)
```
See also: [`Common.rotate_dihedral!`](@ref Common)
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