mutable struct DihedralMutator <: AbstractMutator
    angle_sampler::Function # shoudl return a float in radians
    p_mut::AbstractFloat
    step_size::AbstractFloat
    selection::Opt{AbstractSelection}
end


function (dihedral_mutator::DihedralMutator)(pose::Pose)
    if dihedral_mutator.selection === nothing
        atoms = collect(eachatom(pose.graph))
    else
        atoms = dihedral_mutator.selection(pose, gather = true)
    end
    
    dihedral_mutator(pose, atoms)
end

function (dihedral_mutator::DihedralMutator)(pose::Pose, atoms::Vector{Atom})
    for atom in atoms
        if rand() < dihedral_mutator.p_mut
            ∠ = dihedral_mutator.angle_sampler() * dihedral_mutator.step_size
            pose.state[atom].Δϕ += ∠
            ProtoSyn.request_i2c(pose.state, all = true)
        end
    end
end