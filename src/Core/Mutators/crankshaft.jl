mutable struct CrankshaftMutator <: AbstractMutator
    angle_sampler::Function # shoudl return a float in radians
    p_mut::AbstractFloat
    step_size::AbstractFloat
    selection::Opt{AbstractSelection}
end


function (crankshaft_mutator::CrankshaftMutator)(pose::Pose)
    if crankshaft_mutator.selection === nothing
        atoms = collect(eachatom(pose.graph))
    else
        atoms = crankshaft_mutator.selection(pose, gather = true)
    end
    
    crankshaft_mutator(pose, atoms)
end


function (crankshaft_mutator::CrankshaftMutator)(pose::Pose, atoms::Vector{Atom})
    for (i, atom_i) in enumerate(atoms)
        for atom_j in atoms[(i+1):end]
            if rand() < crankshaft_mutator.p_mut
                # 1) Get angle
                ∠ = crankshaft_mutator.angle_sampler() * crankshaft_mutator.step_size
                println("Angle: $(rad2deg(∠))°")
                # 2) Get axis
                # Note 1: CAi is the pivot
                pivot = pose.state[atom_i].t
                axis  = collect(pose.state[atom_j].t - pivot)
                # 3) Get rotation matrix from axis and angle
                rmat = ProtoSyn.rotation_matrix_from_axis_angle(axis, ∠)
                display(rmat)
                println("\n")
                # 4) Apply the rotation matrix to all atoms between CAi and CAj
                # Note 2: Sidechain of last CA must be included in the rotation
                scj = atom_j.container["C"].id - 1
                M = pose.state.x[:, atom_i.id:scj]
                println("Before: $(pose.state[atom_i.id + 1].t)")
                result = (rmat * (M .- pivot)) .+ pivot
                display(pose.state.x)
                pose.state.x[:, atom_i.id:scj] = result
                println("After: $(pose.state[atom_i.id + 1].t)")
                ProtoSyn.request_c2i(pose.state, all = true)
                return pose.state
            end
        end
    end
end