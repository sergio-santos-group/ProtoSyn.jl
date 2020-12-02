using StatsBase

mutable struct DesignMutator <: AbstractMutator
    p_mut::AbstractFloat
    grammar::LGrammar
    selection::Opt{AbstractSelection}
end


function (design_mutator::DesignMutator)(pose::Pose)
    if design_mutator.selection === nothing
        residues = collect(eachatom(pose.graph))
    else
        sele  = design_mutator.selection
        residues = ProtoSyn.promote(sele, Residue)(pose, gather = true)
    end

    design_mutator(pose, residues)
end

function (design_mutator::DesignMutator)(pose::Pose, residues::Vector{Residue})
    for residue in residues
        if rand() < design_mutator.p_mut

            # 1) Get different aminoacid
            cr_name = Peptides.three_2_one[residue.name]
            nr_name = cr_name
            while nr_name == cr_name
                nr_name = sample(Peptides.available_aminoacids)
            end

            # 2) Perform mutation (already requests i2c)
            derivation = [string(nr_name)]
            Peptides.mutate!(pose, residue, design_mutator.grammar, derivation)
        end
    end
end