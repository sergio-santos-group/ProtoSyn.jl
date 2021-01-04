mutable struct CrankshaftMutator <: AbstractMutator
    angle_sampler::Function # Should return a float in radians
    p_mut::AbstractFloat
    step_size::AbstractFloat
    selection::Opt{AbstractSelection}
    exclude_from_last_res::Opt{AbstractSelection} # an"^C$|^O$"r in peptides
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

    if pose.state.i2c
        sync!(pose)
    end

    for (i, atom_i) in enumerate(atoms)
        for atom_j in atoms[(i+1):end]
            if rand() < crankshaft_mutator.p_mut
                # 1) Get angle
                ∠      = crankshaft_mutator.angle_sampler()
                ∠     *= crankshaft_mutator.step_size

                # 2) Get axis
                # Note 1: Ai is the pivot
                pivot  = pose.state[atom_i].t
                axis   = collect(pose.state[atom_j].t - pivot)

                # 3) Get rotation matrix from axis and angle
                rmat   = ProtoSyn.rotation_matrix_from_axis_angle(axis, ∠)

                # 4) Apply the rotation matrix to all atoms between Ai and Aj
                # Note 2: The whole residue where Aj in included should be
                # included in the rotation. In the case of peptides, the C=O
                # atoms should be excluded. This should be set in the
                # exclude_from_last_res selection.
                scj    = atom_j.container.items[end].id
                lr     = SerialSelection{Residue}(atom_j.container.id, :id)
                elr    = crankshaft_mutator.exclude_from_last_res
                sele   = RangeSelection{Atom}(atom_i:scj, :id) & !(lr & elr)
                idxs   = findall(sele(pose).content)
                M      = pose.state.x[:, idxs]
                result = (rmat * (M .- pivot)) .+ pivot
                pose.state.x[:, atom_i.id:scj] = result
                ProtoSyn.request_c2i(pose.state, all = true)

                return pose.state
            end
        end
    end
end