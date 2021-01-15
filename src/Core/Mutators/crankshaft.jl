using Printf

mutable struct CrankshaftMutator <: AbstractMutator
    angle_sampler::Function # Should return a float in radians
    p_mut::AbstractFloat
    step_size::AbstractFloat
    selection::Opt{AbstractSelection}
    exclude_from_last_res::Opt{AbstractSelection} # an"^C$|^O$"r in peptides
    # Note: the "exclude_from_last_res" parameter is a selection. All atoms in
    # the selection AND in the last residue considered for the crankshaft
    # movement are not included in the rotation.
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
                println("Applying crankshaft")
                # 1) Get angle
                ∠      = crankshaft_mutator.angle_sampler()
                ∠     *= crankshaft_mutator.step_size

                # 2) Get axis
                # Note 1: Ai is the pivot
                # sync!(pose)
                pivot  = pose.state[atom_i].t
                axis   = collect(pose.state[atom_j].t - pivot)

                # 3) Get rotation matrix from axis and angle
                rmat   = ProtoSyn.rotation_matrix_from_axis_angle(axis, ∠)

                # 4) Apply the rotation matrix to all atoms between Ai and Aj
                # Note 2: The whole residue where Aj in included should be
                # included in the rotation. In the case of peptides, the C=O
                # atoms should be excluded. This should be set in the
                # exclude_from_last_res selection.
                ids = ProtoSyn.ids(travel_graph(atom_i, atom_j))
                sidechain_selection = !(an"^CA$|^N$|^C$|^H$|^O$"r)
                atoms = sidechain_selection(atom_j.container, gather = true)
                ids2 = ProtoSyn.ids(atoms)
                idxs = vcat(ids, ids2)
                # scj    = atom_j.container.items[end].id
                # lr     = SerialSelection{Residue}(atom_j.container.id, :id)
                # elr    = crankshaft_mutator.exclude_from_last_res
                # sele   = RangeSelection{Atom}(atom_i.id:scj, :id) & !(lr & elr)
                # idxs   = findall(sele(pose).content)
                M      = pose.state.x[:, idxs]
                result = (rmat * (M .- pivot)) .+ pivot
                pose.state.x[:, idxs] = result
                ProtoSyn.request_c2i(pose.state, all = true)
            end
        end
    end
end

function Base.show(io::IO, cm::CrankshaftMutator)
    println("  Crankshaft:")
    println(io, "+"*repeat("-", 68)*"+")
    @printf(io, "| %-5s | %-25s | %-30s |\n", "Index", "Field", "Value")
    println(io, "+"*repeat("-", 68)*"+")
    @printf(io, "| %-5d | %-25s | %-30s |\n", 1, "angle_sampler", "Function $(cm.angle_sampler)")
    @printf(io, "| %-5d | %-25s | %-30.4f |\n", 2, "p_mut", cm.p_mut)
    @printf(io, "| %-5d | %-25s | %-30.4f |\n", 3, "step_size", cm.step_size)
    if cm.selection === nothing
        @printf(io, "| %-5d | %-25s | %-30s |\n", 4, "selection", "Not set")
    else
        @printf(io, "| %-5d | %-25s | %-30s |\n", 4, "selection", "Set: $(typeof(cm.selection).name)")
    end
    if cm.exclude_from_last_res === nothing
        @printf(io, "| %-5d | %-25s | %-30s |\n", 4, "exclude_from_last_res", "Not set")
    else
        @printf(io, "| %-5d | %-25s | %-30s |\n", 4, "exclude_from_last_res", "Set: $(typeof(cm.selection).name)")
    end
    println(io, "+"*repeat("-", 68)*"+")
end