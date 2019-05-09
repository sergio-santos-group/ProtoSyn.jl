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
    BlockrotMutator(blocks::Vector{Common.BlockMetadata}, angle_sampler::Function, p_mut::Float64, step_size::Float64, n_tries::Int64, loop_closer::Drivers.AbstractDriver)

Holds all necessary parameters for the correct simulation of blockrot movements.

# Arguments
- `blocks::Vector{Common.BlockMetadata}`: List of blocks avaliable to be rotated in crankshaft movements, containing all necessary metadata information.
- `angle_sampler::Function`: Function responsible for defining the rotation angle. Should return a `Float64`.
- `p_mut::Float64`: Probability of rotation of each pair of alpha carbons.
- `step_size::Float64`: Scalar that defines the amount of rotation resulting for a blockrot movement.
- `translation_step_size::Float64`: Scalar that defines the amount of translation resulting for a blockrot movement.
- `n_tries::Int64`: Number of attempts to rotate an accepted block. A rotation is deemed invalid when the block ends rotate to a position where the loop cannot close.
- `loop_closer::Drivers.AbstractDriver`: Driver responsible for closing the loops (e.g. Steepest Descent)

# Examples
```julia-repl
julia> Mutators.Blockrot.BlockrotMutator(metadata.blocks, randn, 0.05, 0.25, 10, loop_closer)
CrankshaftMutator(blocks=6, angle_sampler=randn, p_pmut=0.05, step_size=0.25, n_tries=10)
```
See also: [`run!`](@ref)
"""
mutable struct MutatorConfig{F <: Function}  <: Abstract.MutatorConfig
    blocks::Vector{Common.BlockMetadata}
    angle_sampler::F
    p_mut::Float64
    step_size::Float64
    translation_step_size::Float64
    n_tries::Int64
    loop_closer::Abstract.DriverConfig

    function MutatorConfig(blocks::Vector{Common.BlockMetadata}, angle_sampler::Function, p_mut::Float64, step_size::Float64, translation_step_size::Float64, n_tries::Int64, loop_closer::Abstract.DriverConfig)
        new{Function}(blocks, angle_sampler, p_mut, step_size, translation_step_size, n_tries, loop_closer)
    end
end
Base.show(io::IO, b::MutatorConfig) = print(io, "Blockrot.MutatorConfig(blocks=$(length(b.blocks)), angle_sampler=$(string(b.angle_sampler)), p_mut=$(b.p_mut), step_size=$(b.step_size), translation_step_size=$(b.translation_step_size), n_tries=$(b.n_tries))")


@doc raw"""
    run!(state::Common.State, mutator::CrankshaftMutator)

Iterate over a list of [`Common.BlockMetadata`](@ref) and perform blockrot movements on the current
[`Common.State`](@ref). The probability of each block undergo movement is defined in the
[`BlockrotMutator`](@ref).`p_mut`. The new angle is obtained from [`BlockrotMutator`](@ref).`angle_sampler`, who should
return a `Float64` in radians.
After movement, the [`Common.State`](@ref) is updated with the new conformation.
Returns the number of rotations performed.

# Examples
```julia-repl
julia> Mutators.BlockrotMutator.run!(state, mutator)
```
"""
@inline function apply!(state::Common.State, mutator::MutatorConfig)
    
    count::Int64 = 0
    for (block_index, block) in enumerate(mutator.blocks)
        if rand() < mutator.p_mut
            backup = copy(state.xyz)
            pivot = state.xyz[block.pivot, :]'
            for n_try in 1:mutator.n_tries
                # println("Try $n_try")
                rotation_axis  = Aux.rand_vector_in_sphere()
                # println("   Rotation Axis: $rotation_axis")
                angle = mutator.angle_sampler()
                # println("  Rotation Angle: $angle")
                rmat = Aux.rotation_matrix_from_axis_angle(rotation_axis, angle)
                state.xyz[block.atoms, :] = (rmat * (state.xyz[block.atoms, :] .- pivot)')' .+ pivot # Rotation
                translation_axis  = Aux.rand_vector_in_sphere()
                # println("Translation Axis: $translation_axis")
                state.xyz[block.atoms, :] .+= translation_axis' * mutator.translation_step_size       # Translation
                
                #Check if it's plausible to close
                if block_index > 1
                    d_left = zeros(Float64, 3)
                    d_left = norm(@. state.xyz[mutator.blocks[block_index - 1].connector_right, :] - state.xyz[block.connector_left, :])
                    if d_left > block.range_left
                        # println("Aborted blockrot movement on block $block_index (left): $d_left/$(block.range_left)")
                        state.xyz[:] = backup
                        continue
                    end
                end
                if block_index < length(mutator.blocks)
                    d_right = zeros(Float64, 3)
                    d_right = norm(@. state.xyz[mutator.blocks[block_index + 1].connector_left, :] - state.xyz[block.connector_right, :])
                    if d_right > mutator.blocks[block_index + 1].range_left
                        # println("Aborted blockrot movement on block $block_index (right): $d_right/$(mutator.blocks[block_index + 1].range_left)")
                        state.xyz[:] = backup
                        continue
                    end
                end
                break
            end

            typeof(mutator.loop_closer).name.module.run!(state, mutator.loop_closer)
            # mutator.loop_closer.run!(state, mutator.loop_closer)
            count += 1
        end
    end
    return count
end

end