"""
# TODO Documentation
"""
module GMX

    using ProtoSyn
    using ProtoSyn.Peptides
    using Printf

    """
    # TODO: Documentation
    """
    function generate_gmx_itp(pose::Pose, selection::Opt{AbstractSelection};
        molecule_itp_filename::String  = "protein.itp",
        overwrite::Bool            = true,
        keep_temp_files::Bool      = false,
        gmx_histidine_type::Int    = 0,
        gmx_water_model::String    = "tip3p",
        gmx_forcefield::String     = "amber99sb-ildn")

        ProtoSyn.GMX.check_installation()

        # Generate Pose to file from just the selected region. No check on
        # validity on the given selection is done.
        frag      = Pose(fragment(pose, selection))
        temp_path = tempname(".")
        temp_file = temp_path*".pdb"
        ProtoSyn.write(frag, temp_file)

        # Run acpype to generate topology files
        run(pipeline(`echo "$gmx_histidine_type"`, `gmx pdb2gmx -f "$temp_file" -water $gmx_water_model -ff $gmx_forcefield -his`))

        # Extract the [ moleculetype ] section to a separate .itp file
        # Note that the [ atomtypes ] section is included in the forcefield
        function write_moleculetype_to_file(molecule_itp_filename::String)
            if isfile(molecule_itp_filename)
                if overwrite
                    @warn "File $molecule_itp_filename already exists in the current directory. ProtoSyn will overwrite existing file. Check if this is the desired behaviour. Consider setting the `overwrite` flag to `false` or setting a new `molecule_itp_filename`."
                    open(molecule_itp_filename, "w") do io; end
                else
                    error("File $molecule_itp_filename already exists in the current directory. Consider setting the `overwrite` flag to `true` or setting a new `molecule_itp_filename`.")
                end
            end

            open(molecule_itp_filename, "a") do io_molecule
                open("topol.top", "r") do io
                    in_moleculetype = false
                    to_skip         = 0
                    for line in eachline(io)
                        if to_skip > 0
                            to_skip -= 1
                            continue
                        end
                        if startswith(line, "[ moleculetype ]")
                            in_moleculetype = true
                        elseif startswith(line, "; Include Position restraint file")
                            in_moleculetype = false
                        end

                        if in_moleculetype
                            Base.write(io_molecule, line*"\n")
                            if startswith(line, "; Name")
                                Base.write(io_molecule, @sprintf("%-17s %-d\n", molecule_itp_filename[1:(end-4)], 3))
                                to_skip += 1
                            end
                        end
                    end
                end
            end
        end

        write_moleculetype_to_file(molecule_itp_filename)

        # Remove temporary files
        if !keep_temp_files
            for file in readdir(".")
                endswith(file, "#") && rm(file)
                startswith(file, "jl_") && rm(file, recursive = true)
            end
        end
    end # function


    """
    # TODO: Documentation
    Untested on proteins with multiple ligands.
    """
    function generate_gmx_topology(pose::Pose;
        protein_itp_filename::String   = "protein.itp",
        atomtypes_itp_filename::String = "atomtypes.itp",
        molecule_itp_filename::String  = "molecule.itp",
        overwrite::Bool            = true,
        keep_temp_files::Bool      = false,
        gmx_histidine_type::Int    = 0,
        gmx_water_model::String    = "tip3p",
        gmx_forcefield::String     = "amber99sb-ildn",
        topology_filename::String  = "topol.top",
        include_atomtypes::Bool    = false)
    
        # Attempts to identify protein and non-protein molecules based on
        # ProteinSelection
        # Generate .itp for protein molecule
        generate_gmx_itp(pose, ProteinSelection(),
            molecule_itp_filename  = protein_itp_filename,
            overwrite          = overwrite,
            keep_temp_files    = keep_temp_files,
            gmx_histidine_type = gmx_histidine_type,
            gmx_water_model    = gmx_water_model,
            gmx_forcefield     = gmx_forcefield)

        # Generate .itp for non-protein molecules (if present)
        N_non_protein_molecules = count(!ProteinSelection()(pose))
        if N_non_protein_molecules > 0
            ProtoSyn.GMX.generate_gmx_itp(pose, !ProteinSelection(),
                atomtypes_itp_filename   = atomtypes_itp_filename,
                molecule_itp_filename    = molecule_itp_filename,
                overwrite                = overwrite,
                keep_temp_files          = keep_temp_files)
        end

        # Merge .itps in a single .top
        if isfile(topology_filename)
            if overwrite
                @warn "File $topology_filename already exists in the current directory. ProtoSyn will overwrite existing file. Check if this is the desired behaviour. Consider setting the `overwrite` flag to `false` or setting a new `topology_filename`."
                open(topology_filename, "w") do io; end
            else
                error("File $topology_itp_filename already exists in the current directory. Consider setting the `overwrite` flag to `true` or setting a new `topology_filename`.")
            end
        end
        open(topology_filename, "w") do topol
            Base.write(topol, "; Topology file generated by ProtoSyn\n\n")
            Base.write(topol, "; Include forcefield parameters\n")
            Base.write(topol, "#include \"$(gmx_forcefield).ff/forcefield.itp\"\n")
            if N_non_protein_molecules > 0 && include_atomtypes
                Base.write(topol, "#include \"$(atomtypes_itp_filename)\"\n")
            end
            Base.write(topol, "\n; Include water topology\n")
            Base.write(topol, "#include \"$(gmx_forcefield).ff/$(gmx_water_model).itp\"\n")
            Base.write(topol, "\n; Include topology for ions\n")
            Base.write(topol, "#include \"$(gmx_forcefield).ff/ions.itp\"\n")
            Base.write(topol, "\n; Include topology for protein\n")
            Base.write(topol, "#include \"$(protein_itp_filename)\"\n")
            if N_non_protein_molecules > 0
                Base.write(topol, "\n; Include topology for ligands\n")
                Base.write(topol, "#include \"$(molecule_itp_filename)\"\n")
            end
            Base.write(topol, "\n[ system ]\n; Name\n$(pose.graph.name)\n")
            Base.write(topol, "\n[ molecules ]\n; Compound        #mols\n")
            Base.write(topol, @sprintf("%-19s %-d\n", protein_itp_filename[1:(end-4)], 1))
            if N_non_protein_molecules > 0
                Base.write(topol, @sprintf("%-19s %-d\n", molecule_itp_filename[1:(end-4)], 1))
            end
        end

        if !keep_temp_files
            rm("conf.gro")
            rm("posre.itp")
        end
    end


    struct GMXSA <: ProtoSyn.SearchAlgorithm end

    function (sa::GMXSA)(atom::Atom, stack::Vector{Atom})

        if atom.container.name == "THR" && atom.name == "CB"
            bonds = copy(ProtoSyn.sort_children(atom))
            bonds = bonds[[1, 3, 2]]
        else
            bonds = copy(ProtoSyn.sort_children(atom))
        end

        if atom.name == "CA"
            c_index = findfirst((atom) -> atom.name == "C", bonds)
            if c_index === nothing
                ProtoSyn.verbose.mode && @warn "Tried to sort $atom children following IUPAC convention, but no \"C\" atom was found."
            else
                push!(bonds, bonds[c_index])
                deleteat!(bonds, c_index)
            end
        end
        return vcat(bonds, stack)
    end


    """
    # TODO: Documentation
    """
    function sort_atoms_and_graph_gmx!(pose::Pose)

        # Deal with Prolines
        for container in rn"PRO"(pose, gather = true)

            @assert container.name == "PRO" "Can't sort Proline atoms and graph on a non-Proline residue ($container)."

            # Manually identify necessary heavy-atoms
            N  = container["N"];   N === nothing && error("Tried to sort atoms and graph in residue $container, but no \"N\" atom was found. Check atom names: consider using Peptides.assign_default_atom_names!.")
            CD = container["CD"]; CD === nothing && error("Tried to sort atoms and graph in residue $container, but no \"CD\" atom was found. Check atom names: consider using Peptides.assign_default_atom_names!.")
            CG = container["CG"]; CG === nothing && error("Tried to sort atoms and graph in residue $container, but no \"CG\" atom was found. Check atom names: consider using Peptides.assign_default_atom_names!.")
            CB = container["CB"]; CB === nothing && error("Tried to sort atoms and graph in residue $container, but no \"CB\" atom was found. Check atom names: consider using Peptides.assign_default_atom_names!.")
            CA = container["CA"]; CA === nothing && error("Tried to sort atoms and graph in residue $container, but no \"CA\" atom was found. Check atom names: consider using Peptides.assign_default_atom_names!.")

            # Manually assign parents
            ProtoSyn.popparent!(CD); ProtoSyn.setparent!(CD, N)
            ProtoSyn.popparent!(CG); ProtoSyn.setparent!(CG, CD)
            ProtoSyn.popparent!(CB); ProtoSyn.setparent!(CB, CG)
            ProtoSyn.popparent!(CA); ProtoSyn.setparent!(CA, CB)    
        end
        
        # Sort atoms in the container according to the new parenthoods
        # The GMXSA already deals with threonines
        ProtoSyn.sort_atoms_by_graph!(pose, search_algorithm = GMXSA())

        # Re-calculate internal coordinates based on new parenthoods
        ProtoSyn.request_c2i!(pose.state)
        sync!(pose)

        return pose
    end


    """
    # TODO: Documentation
    """
    function assign_gmx_atom_names!(pose::Pose, selection::Opt{AbstractSelection} = nothing)

        conv = Dict{String, String}(
            "HD11" => "HD1",
            "HD12" => "HD2",
            "HD13" => "HD3",
            "CD1"  => "CD",
        )

        conv_keys = collect(keys(conv))

        if selection === nothing
            sele = TrueSelection{Atom}()
        else
            sele = ProtoSyn.promote(selection, Atom)
        end

        for atom in sele(pose, gather = true)
            if atom.container.name == "ILE" && atom.name in conv_keys
                ProtoSyn.rename!(atom, conv[atom.name])
            end
        end

        ct_sele     = ProtoSyn.promote((sele & CTerminalSelection()), Residue)
        c_terminals = ct_sele(pose, gather = true)
        for c_terminal in c_terminals
            println("Focus on c_terminal: $(c_terminal)")
            if c_terminal["O"] === nothing
                @warn "Tried to rename terminal atoms in residue $c_terminal, but no \"O\" atom was found. Check atom names: consider using Peptides.assign_default_atom_names!.\n        Available atoms: $([a.name for a in c_terminal.items])"
            else
                ProtoSyn.rename!(c_terminal["O"], "OC1")
            end
            if c_terminal["OXT"] === nothing
                @warn "Tried to rename terminal atoms in residue $c_terminal, but no \"OXT\" atom was found. Check atom names: consider using Peptides.assign_default_atom_names!.\n        Available atoms: $([a.name for a in c_terminal.items])"
            else
                ProtoSyn.rename!(c_terminal["OXT"], "OC2")
            end
        end

        return pose
    end


    function generate_gmx_files(pose::Pose;
        protein_itp_filename::String   = "protein.itp",
        output_pdb_filename::String    = "system.pdb",
        atomtypes_itp_filename::String = "atomtypes.itp",
        molecule_itp_filename::String  = "molecule.itp",
        overwrite::Bool            = true,
        keep_temp_files::Bool      = false,
        gmx_histidine_type::Int    = 0,
        gmx_water_model::String    = "tip3p",
        gmx_forcefield::String     = "amber99sb-ildn",
        topology_filename::String  = "topol.top",
        include_atomtypes::Bool    = false)

        p = copy(pose)
        sort_atoms_and_graph_gmx!(p)
        assign_gmx_atom_names!(p)
        generate_gmx_topology(p,
            protein_itp_filename   = protein_itp_filename,
            atomtypes_itp_filename = atomtypes_itp_filename,
            molecule_itp_filename  = molecule_itp_filename,
            overwrite          = overwrite,
            keep_temp_files    = keep_temp_files,
            gmx_histidine_type = gmx_histidine_type,
            gmx_water_model    = gmx_water_model,
            gmx_forcefield     = gmx_forcefield,
            topology_filename  = topology_filename,
            include_atomtypes  = include_atomtypes)

        ProtoSyn.write(p, output_pdb_filename)
    end # function

end # module