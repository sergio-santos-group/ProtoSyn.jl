module Blockrot

using ..Forcefield
using ..Drivers
using ..Common
using ..Aux
using ..Print
using Printf
using LinearAlgebra

mutable struct BlockrotMutator
    blocks::Vector{Common.BlockMetadata}
    angle_sampler::Function
    p_mut::Float64
    step_size::Float64
    loop_closer::Drivers.AbstractDriver

    function BlockrotMutator(blocks::Vector{Common.BlockMetadata}, angle_sampler::Function, p_mut::Float64, step_size::Float64, loop_closer::Drivers.AbstractDriver)
        new(blocks, angle_sampler, p_mut, step_size, loop_closer)
    end
end
Base.show(io::IO, b::BlockrotMutator) = print(io, "BlockrotMutator(blocks=$(b.blocks), angle_sampler=$(string(b.angle_sampler)), p_mut=$(b.p_mut), step_size=$(b.step_size))")


@inline function run!(state::Common.State, mutator::BlockrotMutator)
    
    count::Int64 = 0
    for (block_index, block) in enumerate(mutator.blocks)
        if rand() < mutator.p_mut && block_index > 1
            axis  = Aux.rand_vector_in_sphere()
            angle = mutator.angle_sampler()
            rmat = Aux.rotation_matrix_from_axis_angle(axis, angle)
            pivot = state.xyz[block.pivot, :]'
            state.xyz[block.atoms, :] = (rmat * (state.xyz[block.atoms, :] .- pivot)')' .+ pivot
            
            #Check if it's even possible to close
            if block_index > 1
                d_left = zeros(Float64, 3)
                d_left = norm(@. state.xyz[mutator.blocks[block_index - 1].connector_right, :] - state.xyz[block.connector_left, :])
                if d_left > block.range_left
                    println("$d_left/$(block.range_left)")
                    continue
                end
            end
            if block_index < length(mutator.blocks)
                d_right = zeros(Float64, 3)
                d_right = norm(@. state.xyz[mutator.blocks[block_index + 1].connector_left, :] - state.xyz[block.connector_right, :])
                if d_right > mutator.blocks[block_index + 1].range_left
                    println("$d_right/$(mutator.blocks[block_index + 1].range_left)")
                    continue
                end
            end

            mutator.loop_closer.run!(state, mutator.loop_closer)
            count += 1
        end
    end
    return count
end


end