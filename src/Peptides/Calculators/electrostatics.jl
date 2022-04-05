module Electrostatics
    
    using ProtoSyn
    using ProtoSyn.Peptides

    """
    # TODO
    """
    function assign_default_charges!(pose::Pose, res_lib::LGrammar = Peptides.grammar, selection::Opt{AbstractSelection} = nothing; supress_warn::Bool = false)

        for segment in eachsegment(pose.graph)
            # Check caps
            # N-terminal
            N = ProtoSyn.identify_atom_by_bonding_pattern(segment[1], ["N", "C", "C", "O"])
            if length(N.bonds) > 2 # Expects naming H1, H2 and H3
                pose.state[segment[1, "H1"]].δ = 0.3333
                pose.state[segment[1, "H2"]].δ = 0.3333
                pose.state[segment[1, "H3"]].δ = 0.3333
            end

            last_res_id = segment[end].id
            Cs = ProtoSyn.identify_atom_by_bonding_pattern(segment[end], ["C", "C", "N", "C"])
            if isa(Cs, Atom)
                C = copy(Cs)
            else
                C = [c for c in Cs if length(c.bonds) <= 3][1]
            end
            if length(C.bonds) > 2 
                pose.state[segment[last_res_id, "OXT"]].δ = -0.7800
            end
        end

        return ProtoSyn.Calculators.Electrostatics.assign_default_charges!(
            pose, res_lib, selection, supress_warn = supress_warn)
    end
end