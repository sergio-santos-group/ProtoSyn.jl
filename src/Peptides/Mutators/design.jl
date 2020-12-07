using StatsBase

mutable struct DesignMutator <: AbstractMutator
    p_mut::AbstractFloat
    grammar::LGrammar
    selection::Opt{AbstractSelection}
    searchable_aminoacids::Vector{Char}
end

DesignMutator(p_mut::AbstractFloat, grammar::LGrammar, selection::Opt{AbstractSelection}; searchable_aminoacids::Vector{Char} = deepcopy(Peptides.available_aminoacids)) = begin
    return DesignMutator(p_mut, grammar, selection, searchable_aminoacids)
end

function (design_mutator::DesignMutator)(pose::Pose)
    if design_mutator.selection === nothing
        residues = collect(eachresidue(pose.graph))
    else
        sele  = design_mutator.selection
        residues = ProtoSyn.promote(sele, Residue)(pose, gather = true)
    end

    design_mutator(pose, residues)
    # design_mutator(pose, [pose.graph[1][1]])
end

function (design_mutator::DesignMutator)(pose::Pose, residues::Vector{Residue})
    for residue in residues
        if rand() < design_mutator.p_mut

            # 1) Get different aminoacid
            cr_name = Peptides.three_2_one[residue.name]
            nr_name = cr_name
            while nr_name == cr_name
                nr_name = sample(design_mutator.searchable_aminoacids)
            end

            # 2) Perform mutation (already requests i2c)
            derivation = [string(nr_name)]
            Peptides.mutate!(pose, residue, design_mutator.grammar, derivation)
        end
    end
end