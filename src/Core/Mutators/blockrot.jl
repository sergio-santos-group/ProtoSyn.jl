using ProtoSyn: Driver

@doc """
    BlockRotMutator <: AbstractMutator

The BlockRotMutator uses RotationRigidBodyMutator to rotate 1 or more blocks of
contiguous residues. After rotation, the severed ends of the block are
reconnected using a given `close` function or Driver.

# Fields:

| `axis_sampler`::Function -> Should return a Vector{Float}. This is the axis
over which the rotation takes place (See `ProtoSyn.rand_vector_in_sphere`).

| `angle_sampler`::Function -> Should return a Float. This is the rotation
angle, in radians.

| `pivot_sampler`::Function -> Receives a pose + block selection as input;
should return a Vector{Float}. This is the center from which the axis is applied
for the rotation (See `ProtoSyn.center_of_mass`).

| `p_mut`::AbstractFloat -> The probability of mutation of each given block.

| `step_size`::AbstractFloat -> The magnitude of the rotation (multiplied by the
rotation angle).

| `blocks`::Vector{RangeSelection} -> The list of blocks to be rotated. Each
block must be a RangeSelection{Residue}.

| `close`::Union{Function, Driver} -> The closing function or driver, applied
after rotation. No closure check is performed. Example: `SteepestDescent`

# Example:
```jldoctest
julia> ProtoSyn.Mutators.BlockRotMutator(ProtoSyn.rand_vector_in_sphere, randn, ProtoSyn.center_of_mass, 1.0, 0.2, [rid"27:44"], sd)
...
```
"""
mutable struct BlockRotMutator <: AbstractMutator
    axis_sampler::Function
    angle_sampler::Function
    pivot_sampler::Function
    p_mut::AbstractFloat
    step_size::AbstractFloat
    blocks::Vector{RangeSelection}
    close::Union{Function, Driver}
end


function (blockrot_mutator::BlockRotMutator)(pose::Pose)
    for block in blockrot_mutator.blocks
        @assert ProtoSyn.selection_type(block) == ProtoSyn.Residue "Selection blocks must be of type Residue. Consider using 'promote' function."
        if rand() < blockrot_mutator.p_mut
            # 1. Perform rotation
            rotate = ProtoSyn.Mutators.RotationRigidBodyMutator(
                blockrot_mutator.axis_sampler,
                blockrot_mutator.angle_sampler,
                blockrot_mutator.pivot_sampler,
                blockrot_mutator.step_size,
                block)
            rotate(pose)

            # 2. Close ends
            blockrot_mutator.close(pose)
        end
    end
end