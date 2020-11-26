using Statistics

mutable struct TranslationRigidBodyMutator <: AbstractMutator
    translation_vector_sampler::Function # should return a Vector{Float}
    step_size::AbstractFloat
    selection::Opt{AbstractSelection}
end


function (rigid_body_mutator::TranslationRigidBodyMutator)(pose::Pose)
    if rigid_body_mutator.selection === nothing
        atoms = collect(eachatom(pose.graph))
    else
        sele  = rigid_body_mutator.selection
        atoms = ProtoSyn.promote(sele, Atom)(pose, gather = true)
    end
    
    rigid_body_mutator(pose, atoms)
end


function (rigid_body_mutator::TranslationRigidBodyMutator)(pose::Pose, atoms::Vector{Atom})
    translation_vector   = rigid_body_mutator.translation_vector_sampler()
    translation_vector .*= rigid_body_mutator.step_size
    for atom in atoms
        # Note: Since the t setproperty! is being intercepted, we can't do
        # pose.state[atom].t .+= translation_vector
        pose.state[atom].t = translation_vector .+ pose.state[atom].t
    end
    ProtoSyn.request_c2i(pose.state)
end

# ---
# Rotation Rigid Body

mutable struct RotationRigidBodyMutator{T <: AbstractFloat} <: AbstractMutator
    axis_sampler::Function # should return a Vector{Float}
    angle_sampler::Function # should return a Float
    step_size::T
    rotation_pivot::Opt{Vector{T}} # if nothing, will use COM
    selection::Opt{AbstractSelection}
end


function (rigid_body_mutator::RotationRigidBodyMutator)(pose::Pose)
    axis  = rigid_body_mutator.axis_sampler()
    angle = rigid_body_mutator.angle_sampler() * rigid_body_mutator.step_size
    rmat  = ProtoSyn.rotation_matrix_from_axis_angle(axis, angle)
    mask  = ProtoSyn.promote(rigid_body_mutator.selection, Atom)(pose)
    idxs  = findall(mask.content)
    if rigid_body_mutator.rotation_pivot == nothing
        pivot = mean(pose.state.x[:, idxs], dims = 2)
    else
        pivot = rigid_body_mutator.rotation_pivot
    end
    pose.state.x[:, idxs] = (rmat * (pose.state.x[:, idxs] .- pivot)) .+ pivot

    ProtoSyn.request_c2i(pose.state)
end