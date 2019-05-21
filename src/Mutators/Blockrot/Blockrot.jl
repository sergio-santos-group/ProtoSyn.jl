module Blockrot

using ..Forcefield
using ..Abstract
using ..Drivers
using ..Common
using ..Aux
using ..Print
using Printf
using LinearAlgebra

@doc raw"""
    MutatorConfig(; blocks::Vector{Common.BlockMetadata} = [], angle_sampler::Function = rand, p_mut::Float64 = 0.0, step_size::Float64 = 0.0, translation_step_size::Float64 = 0.0, n_tries::Int64 = 0, loop_closer::Drivers.AbstractDriver = Abstract.NullDriverConfig())

Holds all necessary parameters for the correct simulation of blockrot movements.

# Arguments
- `blocks::Vector{Common.BlockMetadata}`: List of blocks avaliable to be rotated in crankshaft movements, containing all necessary metadata information (Default: empty).
- `angle_sampler::Function`: Function responsible for defining the rotation angle. Should return a `Float64` (Default: rand).
- `p_mut::Float64`: Probability of rotation of each pair of alpha carbons (Default: 0.0).
- `step_size::Float64`: Scalar that defines the amount of rotation resulting for a blockrot movement (Default: 0.0).
- `translation_step_size::Float64`: Scalar that defines the amount of translation resulting for a blockrot movement (Default: 0.0).
- `n_tries::Int64`: Number of attempts to rotate an accepted block. A rotation is deemed invalid when the block ends rotate to a position where the loop cannot close (Default: 0).
- `rot_axis::Union{Symbol, Vector{Float64}}`: Rotation axis around which the rotation of the block will be performed. 
- `trans_axis::Union{Symbol, Vector{Float64}}`: Translation axis by which the translation the block will be performed. 
Both can accept a Vector{Float64} or a Symbol keyword. Available keywords:
:random -> Chose a random vector from an uniform distribution around a sphere;
:longitudinal -> Calulate the longitudinal axis along the block length;
(Default: :random)
- `loop_closer::Drivers.AbstractDriver`: Driver responsible for closing the loops (e.g. [`SteepestDescent`](@ref Drivers))

# Examples
```julia-repl
julia> Mutators.Blockrot.MutatorConfig(blocks = metadata.blocks, angle_sampler = randn, p_mut = 0.05, step_size = 0.25, n_tries = 10, loop_closer)
Blockrot.MutatorConfig(blocks=6, angle_sampler=randn, p_pmut=0.05, step_size=0.25, translation_step_size=0.0, n_tries=10)
```
See also: [`apply!`](@ref)
"""
Base.@kwdef mutable struct MutatorConfig{F <: Function}  <: Abstract.MutatorConfig
    
    # Parameters                               # Defaults
    blocks::Vector{Common.BlockMetadata}       = Vector{Common.BlockMetadata}()
    angle_sampler::F                           = rand
    p_mut::Float64                             = 0.0
    step_size::Float64                         = 0.0
    translation_step_size::Float64             = 0.0
    n_tries::Int64                             = 0
    rot_axis::Union{Symbol, Vector{Float64}}   = :random
    trans_axis::Union{Symbol, Vector{Float64}} = :random
    loop_closer::Abstract.DriverConfig         = Abstract.NullDriverConfig()
end # end struct
Base.show(io::IO, b::MutatorConfig) = print(io, "Blockrot.MutatorConfig(blocks=$(length(b.blocks)), angle_sampler=$(string(b.angle_sampler)), p_mut=$(b.p_mut), step_size=$(b.step_size), translation_step_size=$(b.translation_step_size), n_tries=$(b.n_tries), rot_axis=$(b.rot_axis))")


function get_rotation_axis(state::Common.State, rot_axis::Union{Symbol, Vector{Float64}}, block::Common.BlockMetadata)::Vector{Float64}
    if isa(rot_axis, Vector{Float64})
        return rot_axis
    elseif rot_axis == :random
        return Aux.rand_vector_in_sphere()
    elseif rot_axis == :longitudinal
        rot_axis = [0.0, 0.0, 0.0]
        @inbounds for i in 1:3
            rot_axis[i] += state.xyz[block.connector_right, i] - state.xyz[block.connector_left, i]
        end
        return rot_axis
        # return state.xyz[block.connector_right, :] - state.xyz[block.connector_left, :]
    else
        error("Keyword $(rot_axis) not currently supported by ProtoSyn for Blockrot axis definition.")
    end
end


@doc raw"""
    apply!(state::Common.State, mutator::MutatorConfig)

Iterate over a list of [`Common.BlockMetadata`](@ref) and perform blockrot movements on the current
[`Common.State`](@ref). The probability of each block undergo movement is defined in the
[`MutatorConfig`](@ref).`p_mut`. The new angle is obtained from [`MutatorConfig`](@ref).`angle_sampler`, who should
return a `Float64` in radians.
After movement, the [`Common.State`](@ref) is updated with the new conformation.
Returns the number of rotations performed.

# Examples
```julia-repl
julia> Mutators.BlockrotMutator.apply!(state, mutator)
```
"""
@inline function apply!(state::Common.State, mutator::MutatorConfig)
    
    count::Int64 = 0
    for (block_index, block) in enumerate(mutator.blocks)
        if rand() < mutator.p_mut
            backup = copy(state.xyz)
            pivot = state.xyz[block.pivot, :]'
            for n_try in 1:mutator.n_tries

                rotation_axis = get_rotation_axis(state, mutator.rot_axis, block)

                angle = mutator.angle_sampler()
                rmat = Aux.rotation_matrix_from_axis_angle(rotation_axis, angle)
                state.xyz[block.atoms, :] = (rmat * (state.xyz[block.atoms, :] .- pivot)')' .+ pivot # Rotation
                translation_axis = get_rotation_axis(state, mutator.trans_axis, block)
                state.xyz[block.atoms, :] .+= translation_axis' * mutator.translation_step_size       # Translation
                
                #Check if it's plausible to close
                if block_index > 1
                    d_left = zeros(Float64, 3)
                    d_left = norm(@. state.xyz[mutator.blocks[block_index - 1].connector_right, :] - state.xyz[block.connector_left, :])
                    if d_left > block.range_left
                        state.xyz[:] = backup
                        continue
                    end # end if
                end # end if
                if block_index < length(mutator.blocks)
                    d_right = zeros(Float64, 3)
                    d_right = norm(@. state.xyz[mutator.blocks[block_index + 1].connector_left, :] - state.xyz[block.connector_right, :])
                    if d_right > mutator.blocks[block_index + 1].range_left
                        state.xyz[:] = backup
                        continue
                    end # end if
                end # end if
                break
            end # end for

            Drivers.run!(state, mutator.loop_closer)
            count += 1
        end # end if
    end # end for
    return count
end # end function

end # end module