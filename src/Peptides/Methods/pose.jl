using ProtoSyn.Units

"""
    # TODO
"""
function append_fragment!(pose::Pose{Topology}, residue::Residue, grammar::LGrammar, frag::Pose{Segment}; ss::Opt{SecondaryStructureTemplate} = nothing, op = "α")
    
    ProtoSyn.append_fragment!(pose, residue, grammar, frag, op = op)
    sync!(pose) # Apply changes from `append_fragment!`

    # By default, the grammar operator alpha sets the psi dihedral (between the
    # anchor aminoacid and the first aminoacid of the fragment) to 180°. This
    # usually interfered with the current position of the oxygen atom (in the
    # C=O bond. For that reason, either the oxygen atom or the nitrogen atom
    # need to be rotated. In the current version of this method, the oxygen atom
    # position is rotated to accomodate the fragment appendage, by setting its
    # ϕ value to be 180° rotated from the current psi dihedral.
    current_psi = ProtoSyn.getdihedral(pose.state, residue.children[1]["N"])
    target_dihedral = current_psi - 180°
    pose.state[residue["O"]].ϕ = pose.state[residue["C"]].Δϕ - target_dihedral
    ProtoSyn.request_i2c!(pose.state, all = true)
    
    if ss !== nothing
        residues = residue.container.items[(residue.index + 1):(residue.index + length(frag.graph))]
        ProtoSyn.Peptides.setss!(pose.state, ss, residues)
    end

    return pose
end


"""
    # TODO
"""
function append_fragment!(pose::Pose{Topology}, residue::Residue, grammar::LGrammar, derivation; ss::Opt{SecondaryStructureTemplate} = nothing, op = "α")
    frag = ProtoSyn.fragment(grammar, derivation)
    return ProtoSyn.Peptides.append_fragment!(pose, residue, grammar, frag, ss = ss,
        op = op)
end


"""
    # TODO
"""
function insert_fragment!(pose::Pose{Topology}, residue::Residue,
    grammar::LGrammar, frag::Pose{Segment};
    ss::Opt{SecondaryStructureTemplate} = nothing, op = "α")

    # Setup
    frag_size = length(eachresidue(frag.graph))
    anchor    = residue.parent

    # Case being connected to root, perform soft uncap the N-terminal
    anchor == ProtoSyn.root(pose.graph).container && begin
        residue["H2"] !== nothing && ProtoSyn.pop_atom!(pose, residue["H2"])
        residue["H3"] !== nothing && ProtoSyn.pop_atom!(pose, residue["H3"])

        # Rename H1 to H
        if residue["H1"] !== nothing
            residue.itemsbyname["H"] = residue["H1"]
            pop!(residue.itemsbyname, "H1")
            residue["H"].name = "H"
        end
    end

    # Insert the fragments
    ProtoSyn.insert_fragment!(pose, residue, grammar, frag, op = op)

    if !(anchor == ProtoSyn.root(pose.graph).container)
        # If not connected to root, adjust the position of the upstream C=O bond
        current_psi = ProtoSyn.getdihedral(pose.state, anchor.children[1]["N"])
        target_dihedral = current_psi - 180°
        pose.state[anchor["O"]].ϕ = pose.state[anchor["C"]].Δϕ - target_dihedral
        ProtoSyn.request_i2c!(pose.state, all = true)
    else
        # If connected to root, adjust the position of the downstream N=H bond
        current_omega = ProtoSyn.getdihedral(pose.state, anchor.children[1]["CA"])
        target_dihedral = current_omega - 180°
        pose.state[residue["H"]].ϕ = pose.state[anchor.children[1]["N"]].Δϕ - target_dihedral
        ProtoSyn.request_i2c!(pose.state, all = true)
    end

    # Case being connected to root, perform soft uncap the N-terminal
    anchor == ProtoSyn.root(pose.graph).container && begin

    end

    # Apply secondary structure
    if ss !== nothing
        r_id = residue.index
        residues = residue.container.items[r_id:r_id + frag_size]
        ProtoSyn.Peptides.setss!(pose.state, ss, residues)
    end
end


"""
    # TODO
"""
function insert_fragment!(pose::Pose{Topology}, residue::Residue, grammar::LGrammar, derivation; ss::Opt{SecondaryStructureTemplate} = nothing, op = "α")
    frag = ProtoSyn.fragment(grammar, derivation)
    return ProtoSyn.Peptides.insert_fragment!(pose, residue, grammar, frag, ss = ss,
        op = op)
end