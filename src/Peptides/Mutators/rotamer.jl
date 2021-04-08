using ProtoSyn.Peptides.Rotamers
using ProtoSyn.Peptides: Dihedral
using ProtoSyn: TrueSelection, getdihedral

# TODO Documentation

mutable struct RotamerMutator <: AbstractMutator
    rotamer_library::Dict{String, ProtoSyn.Peptides.Rotamers.BBD_RotamerLibrary}
    p_mut::AbstractFloat
    n_first::Int
    selection::Opt{AbstractSelection}
end


function (rotamer_mutator::RotamerMutator)(pose::Pose)
    if rotamer_mutator.selection === nothing
        atoms = collect(eachatom(pose.graph))
    else
        sele = rotamer_mutator.selection
        atoms = ProtoSyn.promote(sele, Atom)(pose, gather = true)
    end

    rotamer_mutator(pose, atoms)
end


function (rotamer_mutator::RotamerMutator)(pose::Pose, atoms::Vector{Atom})
    for atom in atoms
        if rand() < rotamer_mutator.p_mut

            # 1) Get the residue name
            residue = atom.container
            name    = residue.name

            # 2) Get the residue phi & psi
            phi = getdihedral(pose.state, Dihedral.phi(residue))
            psi = getdihedral(pose.state, Dihedral.psi(residue))

            # 3) Sample with correct name, phi & psi
            rotamer_stack = nothing
            try
                rotamer_stack = rotamer_mutator.rotamer_library[name][phi, psi]
            catch KeyError
                continue
            end
            rotamer = Rotamers.sample(rotamer_stack, rotamer_mutator.n_first)

            # 4) Apply sampled rotamer
            Rotamers.apply!(pose.state, rotamer, residue)
            ProtoSyn.request_i2c!(pose.state, all = true)
        end
    end
end