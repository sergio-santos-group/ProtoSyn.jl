module Sidechain

using ..Aux
using ..Common
using ..Abstract
using StatsBase
using Random

@doc raw"""
    MutatorConfig(; sidechains::Vector{Common.Dihedral} = [], rot_lib::Dict{String, Any} = Dict{String, Any}(), p_mut::Float64 = 0.0) <: Abstract.MutatorConfig

Holds all necessary parameters for the correct simulation of sidechain movements.

# Arguments
- `sidechains::Vector{Vector{Common.Dihedral}}`: List of sidechains avaliable to be mutated in sidechain movements (Default = empty).
- `rot_lib::Dict{String, Any}`: Dictionary containing all possible dihedrals for each chi of an aminoacid (Default = empty). 
- `p_mut::Float64`: Probability of mutation of each sidechain (Default = 0.0).

# Examples
```julia-repl
julia> Mutators.Sidechain.MutatorConfig(sidechains = sidechains, rot_lib = rot_lib, p_mut = 0.25)
Sidechain.MutatorConfig(sidechains=68, (...), p_pmut=0.05)
```
See also: [`apply!`](@ref)
"""
Base.@kwdef mutable struct MutatorConfig <: Abstract.MutatorConfig
    
    # Parameters                                # Defaults
    sidechains::Vector{Vector{Common.Dihedral}} = Vector{Vector{Common.Dihedral}}()
    rot_lib::Dict{String, Any}                  = Dict{String, Any}()
    p_mut::Float64                              = 0.0
end # end struct
Base.show(io::IO, b::MutatorConfig) = print(io, "Sidechain.MutatorConfig(sidechains=$(length(b.sidechains)), rot_lib=$(b.rot_lib), p_mut=$(b.p_mut)")


@doc raw"""
    apply!(state::Common.State, mutator::MutatorConfig)

Iterate over a list of [`Dihedral`](@ref Common) and perform sidechain mutations on the current
[`Common.State`](@ref). The probability of each sidechain undergo movements is defined in the
[`MutatorConfig`](@ref).`p_mut`. The new list of angles (one of each `chi` in the sidechain) is obtained
from `rot_lib`, and sampled with a random weighted probability, defined in
`rot_lib` aswell. After mutation, the [`Common.State`](@ref) is updated with the new conformation.
Returns the number of mutations performed.

# Examples
```julia-repl
julia> Mutators.Sidechain.apply!(state, mutator)
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
            end # end for
            count += 1
        end # end if
    end # end for
    return count
end # end function

end # end module