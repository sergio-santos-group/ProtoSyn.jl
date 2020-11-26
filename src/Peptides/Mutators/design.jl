mutable struct DesignMutator <: AbstractMutator
    p_mut::AbstractFloat
    selection::Opt{AbstractSelection}
end


function (design_mutator::DesignMutator)(pose::Pose)
    if design_mutator.selection === nothing
        mask = TrueSelection{Atom}()(pose)
    else
        mask = design_mutator.selection(pose)
    end

    design_mutator(pose, mask)
end

function (design_mutator::RotamerMutator)(pose::Pose, mask::ProtoSyn.Mask{Atom})
    for (index, atom) in enumerate(eachatom(pose.graph))
        if mask[index] && rand() < design_mutator.p_mut

            # 1) Get the residue name
            residue = atom.container
            name = residue.name

            # 2) Get the residue phi & psi
            phi = getdihedral(pose.state, Dihedral.phi(residue))
            psi = getdihedral(pose.state, Dihedral.phi(residue))

            # 3) Sample with correct name, phi & psi
            design_stack = nothing
            try
                design_stack = design_mutator.design_library[name][phi, psi]
            catch KeyError
                continue
            end
            rotamer = Rotamers.sample(design_stack, design_mutator.n_first)

            # 4) Apply sampled rotamer
            Rotamers.apply!(pose.state, rotamer, residue)
            ProtoSyn.request_i2c(pose.state, all = true)
        end
    end
end