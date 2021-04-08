mutable struct TranslationRigidBodyMutator <: AbstractMutator
    translation_vector_sampler::Function # should return a Vector{Float}
    step_size::AbstractFloat
    selection::Opt{AbstractSelection}
end

TranslationRigidBodyMutator(translation_vector_sampler::Function, step_size::AbstractFloat) = begin
    TranslationRigidBodyMutator(translation_vector_sampler, step_size, nothing)
end

function (rigid_body_mutator::TranslationRigidBodyMutator)(pose::Pose)
    if rigid_body_mutator.selection === nothing
        atoms = Vector{Atom}(collect(eachatom(pose.graph)))
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
    ProtoSyn.request_c2i!(pose.state)
end

# ------------------------------------------------------------------------------
# Rotation Rigid Body

mutable struct RotationRigidBodyMutator{T <: AbstractFloat} <: AbstractMutator
    axis_sampler::Function # should return a Vector{Float}
    angle_sampler::Function # should return a Float
    pivot_sampler::Function # receives a pose and a vector of selected
    # atom indexes and should return a Vector{Float}.
    step_size::T
    selection::Opt{AbstractSelection}
end

RotationRigidBodyMutator(axis_sampler::Function, angle_sampler::Function, pivot_sampler::Function, step_size::T) where {T <: AbstractFloat} = begin
    RotationRigidBodyMutator(axis_sampler, angle_sampler, pivot_sampler, step_size, nothing)
end

function (rigid_body_mutator::RotationRigidBodyMutator)(pose::Pose)
    axis  = rigid_body_mutator.axis_sampler()
    angle = rigid_body_mutator.angle_sampler() * rigid_body_mutator.step_size
    rmat  = ProtoSyn.rotation_matrix_from_axis_angle(axis, angle)
    selec = rigid_body_mutator.selection === nothing ? TrueSelection{Atom}() : rigid_body_mutator.selection
    mask  = ProtoSyn.promote(selec, Atom)(pose)
    idxs  = findall(mask.content)
    pivot = rigid_body_mutator.pivot_sampler(pose, idxs)
    pose.state.x[:, idxs] = (rmat * (pose.state.x[:, idxs] .- pivot)) .+ pivot
    ProtoSyn.request_c2i!(pose.state)
end