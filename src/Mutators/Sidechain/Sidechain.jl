module Sidechain

using ..Aux
using ..Common
using ..Abstract
using StatsBase
using Random

@doc raw"""
    SidechainMutator(sidechains::Vector{Common.SidechainMetadata}, rot_lib::Dict{String, Any}, p_mut::Float64)

Holds all necessary parameters for the correct simulation of sidechain movements.

# Arguments
- `sidechains::Vector{Common.SidechainMetadata}`: List of sidechains avaliable to be mutated in sidechain movements.
- `rot_lib::Dict{String, Any}`: Dictionary containing all possible dihedrals for each chi of an aminoacid. 
- `p_mut::Float64`: Probability of mutation of each sidechain.

# Examples
```julia-repl
julia> Mutators.Sidechain.SidechainMutator(sidechains, rot_lib 0.25)
SidechainMutator(sidechains=68, (...), p_pmut=0.05)
```
See also: [`run!`](@ref)
"""
@Base.kwdef mutable struct MutatorConfig <: Abstract.MutatorConfig
    sidechains::Vector{Vector{Common.Dihedral}} = Vector{Vector{Common.Dihedral}}()
    rot_lib::Dict{String, Any} = Dict{String, Any}()
    p_mut::Float64 = 0.0
end
Base.show(io::IO, b::MutatorConfig) = print(io, "Sidechain.MutatorConfig(sidechains=$(length(b.sidechains)), rot_lib=$(b.rot_lib), p_mut=$(b.p_mut)")


@doc raw"""
    run!(state::Common.State, mutator::SidechainMutator)

Iterate over a list of [`SidechainMetadata`](@ref Common.SidechainMetadata) and perform sidechain mutations on the current
[`Common.State`](@ref). The probability of each sidechain undergo movements is defined in the
[`SidechainMutator`](@ref).`p_mut`. The new list of angles (one of each `chi` in the sidechain) is obtained
from [`SidechainMetadata`](@ref Common.SidechainMetadata).`rotamers`, and sampled with a random weighted probability, defined in
[`SidechainMetadata`](@ref Common.SidechainMetadata).`weights`. After mutation, the [`Common.State`](@ref) is updated with the new conformation.
Returns the number of mutations performed.

# Examples
```julia-repl
julia> Mutators.Sidechain.run!(state, mutator)
```
See also: [`rotate_dihedral_to!`](@ref Common.rotate_dihedral_to!)
"""
@inline function apply!(state::Common.State, mutator::MutatorConfig)
    
    count::Int64 = 0
    for sidechain in mutator.sidechains
        if rand() < mutator.p_mut
            aa = Aux.conv321(sidechain[1].residue.name)
            rotamer_index = sample(1:length(mutator.rot_lib[aa]["w"]), Weights(convert(Vector{Float64}, mutator.rot_lib[aa]["w"])))
            rotamer = mutator.rot_lib[aa]["rot"][rotamer_index]
            for chi_index in 1:length(rotamer)
                chi_center = rotamer[chi_index]
                chi_range  = mutator.rot_lib[aa]["range"][rotamer_index][chi_index]
                chi_value  = (chi_range/1.175) * randn() + chi_center
                Common.rotate_dihedral_to!(state.xyz, sidechain[chi_index], chi_value)
            end
            count += 1
        end
    end
    return count
end

end