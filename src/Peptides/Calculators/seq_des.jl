module SeqDes

    using ProtoSyn

    atom_types = Dict{String, Int}(
        "N" => 1,
        "C" => 2,
        "O" => 3,
        "S" => 4,
        "P" => 5,
    )

    residue_types = Dict{String, Int}(
        "HIS" => 0,
        "LYS" => 1,
        "ARG" => 2,
        "ASP" => 3,
        "GLU" => 4,
        "SER" => 5,
        "THR" => 6,
        "ASN" => 7,
        "GLN" => 8,
        "ALA" => 9,
        "VAL" => 10,
        "LEU" => 11,
        "ILE" => 12,
        "MET" => 13,
        "PHE" => 14,
        "TYR" => 15,
        "TRP" => 16,
        "PRO" => 17,
        "GLY" => 18,
        "CYS" => 19,
        "MSE" => 13,
    )

    """
    # TODO: DOCUMENTATION
    Should return:
        - atom_coords: N x 3 Vector (atomic positions)
        - atom_data: N x 4 Vector(atomic data)
            • [residue_id, bb_ind, atom_type, residue_type], where `residue_id`
            is 0-indexed; `bb_ind` is either 1 or 0, if the corresponding atom
            is in the protein backbone or not, respectively; `atom_type` is the
            index of the atom element in the `Calculators.SeqDes.atom_types`
            dictionary; and `residue_type` is the index of the residue type in
            the  `Calculators.SeqDes.residue_types` dictionary.
        - residue_bb_index: Nr x 4 (atomic index of all backbone atoms + CB)
            • [N, CA, C, CB], where all atomic indexes are 0-indexed; for
            residues without CB, the -1 index is used.
        - residue_data: Dict('chain_code' => Nr x 4) (residue data)
            • [residue_id, residue_id_code, residue_index, residue_type], where
            `residue_id_code` can be ignored and set to ' ', residue_index is
            0-indexed and `residue_type` is the index of the residue type in the
            `Calculators.SeqDes.residue_types` dictionary.

        Note: `N` is the number of atoms and `Nr` is the number of residues.
    """
    function get_pdb_data(pose::Pose, selection::Opt{AbstractSelection} = nothing)

        function get_residue_type(residue_name::ProtoSyn.ResidueName)
            if residue_name in keys(residue_types)
                return residue_types[residue_name]
            else
                return 20
            end
        end

        sidechain = SidechainSelection()(pose)

        if selection !== nothing
            atom_sele = ProtoSyn.promote(selection, Atom) & !as"H"
        else
            atom_sele = !as"H"
        end

        
        atoms       = atom_sele(pose, gather = true)
        residues    = collect(eachresidue(pose.graph))
        
        # Atomic descriptors (atom_coords & atom_data)
        T           = eltype(pose.state)
        atom_coords = Vector{Vector{Vector{T}}}()
        atom_data   = Vector{Vector{Vector{Union{Int}}}}()

        # Residue descriptors (residue_bb_index)
        _residue_bb_index = Dict{Int, Dict{String, Int}}()

        for residue in residues
            _residue_bb_index[residue.id] = Dict{String, Int}(
                "N" => -1, "CA" => -1, "C" => -1, "CB" => -1)
        end

        # Atom loop
        for (i, atom) in enumerate(atoms)
            # Note: i is the non-hydrogen index of the atom

            # Atom_coords
            push!(atom_coords, [Vector{T}(pose.state[atom].t)])

            # Atom_data
            residue_id       =  atom.container.id - 1
            bb_ind           = !sidechain[atom.id]
            if atom.symbol in keys(atom_types)
                atom_type    = atom_types[atom.symbol]
            else
                atom_type    = 6
            end
            residue_type = get_residue_type(atom.container.name)
            push!(atom_data, [Int[residue_id, bb_ind, atom_type, residue_type]])

            # Residue_bb_index (1/2)
            if atom.name in ["N", "CA", "C", "CB"]
                _residue_bb_index[atom.container.id][atom.name] = i - 1
            end

        end

        # Residue loop
        residue_bb_index = Vector{Vector{Int}}()
        residue_data = Dict{String, Vector{Tuple{Int, String, Int, Int}}}()
        for segment in eachsegment(pose.graph)
            nv = Vector{Tuple{Int, String, Int, Int}}()
            residue_data[string(segment.code)] = nv
        end

        for residue in residues
            
            # Residue_bb_index (2/2)
            r = _residue_bb_index[residue.id]
            push!(residue_bb_index, Int[r["N"], r["CA"], r["C"], r["CB"]])

            # Residue_data
            residue_id      = residue.id
            residue_id_code = " " # ?
            residue_index   = residue.id - 1
            residue_type    = get_residue_type(residue.name)
            push!(residue_data[string(residue.container.code)],
                (residue_id, residue_id_code, residue_index, residue_type))
        end

        return atom_coords, atom_data, residue_bb_index, residue_data
    end

end