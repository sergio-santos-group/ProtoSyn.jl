module Electrostatics
    
    using ProtoSyn
    using ProtoSyn.Peptides

    """
        assign_default_charges!(pose::Pose, [res_lib::LGrammar = Peptides.grammar], [selection::Opt{AbstractSelection}]; [supress_warn::Bool = false])

    Assign default charges to [`Pose`](@ref) `pose` from the given
    [`LGrammar`](@ref) `res_lib` entry, by [`Atom`](@ref) name. If an
    `AbstractSelection` `selection` is provided, only apply charges to the
    selected [`Atom`](@ref) instances. For non-canonical aminoacids and ligands
    (any [`Residue`](@ref) without an entry on `ProtoSyn.three_2_one`
    dictionary) and any [`Residue`](@ref) whose template have different
    [`Atom`](@ref) names, a warning is shown. Set `supress_warn` to `true` to
    ignore these warnings (`false`, by default).

    !!! ukw "Note:"
        This method is an overload of [`assign_default_charges!`](@ref ProtoSyn.Calculators.Electrostatics.assign_default_charges!). Check it's documentation page for more information.
        In comparison, this method is specialized for Peptides (automatically identifies terminal atoms and adds correct charges their caps).

    # See also
    [`assign_acc2_eem_charges_from_file!`](@ref ProtoSyn.Calculators.Electrostatics.assign_acc2_eem_charges_from_file!)

    # Examples
    ```jldoctest
    julia> ProtoSyn.Peptides.Calculators.Electrostatics.assign_default_charges!(pose)
    1143-element Vector{Float64}:
     -0.025115728872692304
     -0.025115728872692304
     -0.025115728872692304
     (...)
    ```
    """
    function assign_default_charges!(pose::Pose, res_lib::LGrammar = Peptides.grammar, selection::Opt{AbstractSelection} = nothing; supress_warn::Bool = false)

        if selection !== nothing
            segments = ProtoSyn.promote(selection, Segment)(pose, gather = true)
        else
            segments = collect(eachsegment(pose.graph))
        end
        
        for segment in segments
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
            if length(C.bonds) > 2 && segment[last_res_id, "OXT"] !== nothing
                pose.state[segment[last_res_id, "OXT"]].δ = -0.7800
            end
        end

        return ProtoSyn.Calculators.Electrostatics.assign_default_charges!(
            pose, res_lib, selection, supress_warn = an"H1|H2|H3|OXT"r)
    end
end