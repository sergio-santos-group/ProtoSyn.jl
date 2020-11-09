mutable struct DihedralMutator <: AbstractMutator
    angle_sampler::Function
    p_mut::AbstractFloat
    step_size::AbstractFloat
    selection::Opt{AbstractSelection}
end


function (dihedral_mutator::DihedralMutator)(pose::Pose)
    if dihedral_mutator.selection === nothing
        mask = TrueSelection{Atom}()(pose)
    else
        mask = dihedral_mutator.selection(pose)
    end
    
    dihedral_mutator(pose, mask)
end


function (dihedral_mutator::DihedralMutator)(pose::Pose, mask::ProtoSyn.Mask{Atom})
    for (index, atom) in enumerate(eachatom(pose.graph))
        if mask[index] && rand() < dihedral_mutator.p_mut
            ∠ = dihedral_mutator.angle_sampler() * dihedral_mutator.step_size
            pose.state[atom].Δϕ += ∠
            ProtoSyn.request_i2c(pose.state, all = true)
        end
    end
end