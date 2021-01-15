using Printf

mutable struct DihedralMutator <: AbstractMutator
    angle_sampler::Function # Should return a float in radians
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

function Base.show(io::IO, dm::DihedralMutator)
    println("  Dihedral:")
    println(io, "+"*repeat("-", 68)*"+")
    @printf(io, "| %-5s | %-25s | %-30s |\n", "Index", "Field", "Value")
    println(io, "+"*repeat("-", 68)*"+")
    @printf(io, "| %-5d | %-25s | %-30s |\n", 1, "angle_sampler", "Function $(dm.angle_sampler)")
    @printf(io, "| %-5d | %-25s | %-30.4f |\n", 2, "p_mut", dm.p_mut)
    @printf(io, "| %-5d | %-25s | %-30.4f |\n", 3, "step_size", dm.step_size)
    if dm.selection === nothing
        @printf(io, "| %-5d | %-25s | %-30s |\n", 4, "selection", "Not set")
    else
        @printf(io, "| %-5d | %-25s | %-30s |\n", 4, "selection", "Set: $(typeof(dm.selection).name)")
    end
    println(io, "+"*repeat("-", 68)*"+")
end