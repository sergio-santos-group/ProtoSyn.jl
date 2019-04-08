module Sidechain

using StatsBase
using Random
using ..Common
using ..Aux

rng = MersenneTwister(1234)

# @doc raw"""
#     DihedralMutator(dihedrals::Vector{Common.Dihedral}, angle_sampler::Function, p_mut::Float64, step_size::Float64)

# Holds all necessary parameters for the correct simulation of dihedral movements.

# # Arguments
# - `dihedrals::Vector{Common.Dihedral}`: List of dihedrals avaliable to be rotated in dihedral movements.
# - `angle_sampler::Function`: Function responsible for defining the new angle for the dihedral. Should return a `Float64`.
# - `p_mut::Float64`: Probability of rotation of each dihedral.
# - `step_size::Float64`: Scalar that defines the amount of change resulting for a dihedral movement.

# # Examples
# ```julia-repl
# julia> Mutators.Diehdral.DihedralMutator(dihedrals, 0.05, randn, 0.25)
# DihedralMutator(dihedrals=68, p_pmut=0.05, angle_sampler=randn, step_size=0.25)
# ```
# See also: [`run!`](@ref)
# """
mutable struct SidechainMutator
    sidechains::Vector{Vector{Common.Dihedral}}
    rot_lib::Dict{String, Dict{String, Any}}
    p_mut::Float64
end
# SidechainMutator(dihedrals::Vector{Common.Dihedral}, angle_sampler::Function; p_mut = 0.0, step_size = 0.0) = DihedralMutator(dihedrals, angle_sampler, p_mut, step_size)
# Base.show(io::IO, b::DihedralMutator) = print(io, "DihedralMutator(dihedrals=$(length(b.dihedrals)), angle_sampler=$(string(b.angle_sampler)), p_mut=$(b.p_mut), step_size=$(b.step_size))")


# @doc raw"""
#     run!(state::Common.State, mutator::DihedralMutator)

# Iterate over a list of [`Common.Dihedral`](@ref) (`dihedrals`) and perform dihedral movements on the current
# [`Common.State`](@ref). The probability of each dihedral undergo movements is defined in the
# [`DihedralMutator`](@ref).`p_mut`. The new angle is obtained from [`DihedralMutator`](@ref).`angle_sampler`, who should
# return a `Float64` in radians.
# After movement, the [`Common.State`](@ref) is updated with the new conformation.
# Returns the number of rotations performed.

# # Examples
# ```julia-repl
# julia> Mutators.Dihedral.run!(state, mutator)
# ```
# See also: [`Common.rotate_dihedral!`](@ref Common)
# """
@inline function run!(state::Common.State, mutator::SidechainMutator)
    
    count::Int64 = 0
    for sidechain in mutator.sidechains
        if rand() < mutator.p_mut
            aa = Aux.conv321(sidechain[1].residue.name)
            weights = mutator.rot_lib[aa]["w"]
            rotamer_index = sample(1:length(weights), Weights(convert(Array{Float64, 1}, weights)))
            for chi_index in 1:mutator.rot_lib[aa]["n_chis"]
                chi_center = mutator.rot_lib[aa]["rot"][rotamer_index][chi_index]
                chi_range  = mutator.rot_lib[aa]["range"][rotamer_index][chi_index]
                chi_value  = (chi_range/1.175) * randn(rng, Float64, 1) .+ chi_center
                Common.rotate_dihedral_to!(state.xyz, sidechain[chi_index], chi_value[1])
            end
            count += 1
        end
    end
    return count
end

end