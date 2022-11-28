"""
    GMX (module)

Holds auxiliary methods for generating .itp, .top and other necessary files for
GROMACS simulations.
"""
module GMX

    using ProtoSyn
    using ProtoSyn.Peptides
    using Printf

    """
        generate_gmx_itp(pose::Pose, [selection::Opt{AbstractSelection} = nothing]; [molecule_itp_filename::String = "protein.itp"], [overwrite::Bool = true], [keep_temp_files::Bool = false], [gmx_histidine_type::Int = 0], [gmx_water_model::String = "tip3p"], [gmx_forcefield::String = "amber99sb-ildn"])

    Generate a protein .itp file for usage in Gromacs using the `gmx pdb2gmx`
    utility program with forcefield parameters for the given [`Pose`](@ref)
    `pose`. If an `AbstractSelection` `selection` is provided, only the selected
    [`Atom`](@ref) instances will be considered (in a [`Fragment`](@ref)) for
    .itp generation. Note that since `gmx pdb2gmx` is being used, non-protein
    [`Residue`](@ref) instance may cause problems or undesired parameters in the
    resulting .itp file. It is reccomended to select protein [`Residue`](@ref)
    instances only. The .itp data will be written to `molecule_itp_filename` (if
    file exists, setting `overwrite` to `true` will overwrite any existing file
    with the same name, set to `true` by default). Setting `keep_temp_files` to
    `true` (`false`, by default) prevents ProtoSyn from deleting any temporary
    files and directories generated during the method call (this includes
    Gromacs backups ending with "#" and any file or directory starting with
    "jl_"). During `gmx pdb2gmx` call, automatically sets:

    * Histidines to be of type `gmx_histidine_type` (default: 0);
    * Water model to be of type `gmx_water_model` (default: "tip3p");
    * Forcefield to be of type `gmx_forcefield` (default: "amber99sb-ildn").

    # See also
    [`generate_gmx_topology`](@ref)

    # Examples
    ```
    julia> ProtoSyn.Peptides.GMX.generate_gmx_itp(pose)
    ✓ All necessary packages were found!
     (...)
    ```
    """
    function generate_gmx_itp(pose::Pose, selection::Opt{AbstractSelection} = nothing;
        molecule_itp_filename::String  = "protein.itp",
        overwrite::Bool            = true,
        keep_temp_files::Bool      = false,
        gmx_histidine_type::Int    = 0,
        gmx_water_model::String    = "tip3p",
        gmx_forcefield::String     = "amber99sb-ildn")

        ProtoSyn.GMX.check_installation()

        if selection === nothing
            selection = TrueSelection{Atom}()
        end

        # Generate Pose to file from just the selected region. No check on
        # validity on the given selection is done.
        if count(selection(pose)) === 0
            @warn "The provided selection yielded no Residue instances for .itp generation using PDB2GMX."
            return nothing
        end

        frag      = Pose(fragment(pose, selection))
        temp_path = tempname(".")
        temp_file = temp_path*".pdb"
        ProtoSyn.write(frag, temp_file)

        # Run acpype to generate topology files
        N_histidines = count(rn"HIS"(pose))
        gmx_histidine_cmd = Base.join(["$gmx_histidine_type" for i in 1:N_histidines], " \n ")
        run(pipeline(`echo "$gmx_histidine_cmd"`, `gmx pdb2gmx -f "$temp_file" -water $gmx_water_model -ff $gmx_forcefield -his`))

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
        generate_gmx_topology(pose::Pose; [protein_itp_filename::String = "protein.itp"], [atomtypes_itp_filename::String = "atomtypes.itp"], [molecule_itp_filename::String = "molecule.itp"], [overwrite::Bool = true], [keep_temp_files::Bool = false], [gmx_histidine_type::Int = 0], [gmx_water_model::String = "tip3p"], [gmx_forcefield::String = "amber99sb-ildn"], [topology_filename::String = "topol.top"], [include_atomtypes::Bool = false])
    
    Generate a system .top file for usage in Gromacs (using the `gmx pdb2gmx`
    and `acpype` utility programs) with forcefield parameters for the given
    [`Pose`](@ref) `pose`. This method automatically tries to select protein
    [`Residue`](@ref) instances using the [`ProteinSelection`](@ref). The
    selected [`Residue`](@ref) instances forcefield parameters will be compiled
    into `protein_itp_filename` using the `gmx pdb2gmx` utility program (if
    file exists, setting `overwrite` to `true` will overwrite any existing file
    with the same name, set to `true` by default). During `gmx pdb2gmx` call,
    automatically sets:

    * Histidines to be of type `gmx_histidine_type` (default: 0);
    * Water model to be of type `gmx_water_model` (default: "tip3p");
    * Forcefield to be of type `gmx_forcefield` (default: "amber99sb-ildn").

    Any remaining [`Residue`](@ref) instances will be considered as ligand
    molecules, and the forcefield information will be compiled into
    `molecule_itp_filename` and `atomtypes_itp_filename` using the `acpype`
    utility program. Setting `keep_temp_files` to `true` (`false`, by default)
    prevents ProtoSyn from deleting any temporary files and directories
    generated during the method call (this includes Gromacs backups ending with
    "#" and any file or directory starting with "jl_"). The resulting .itp files
    are combined into the final `topology_filename` .top file (if
    `include_atomtypes` is set to `true`, the ligand molecule's atomtypes file
    is also included in the topology, `false` by default).

    !!! ukw "Note:"
        This method is untested on proteins with multiple ligands and may fail. This may change in future versions of ProtoSyn.

    # See also
    [`generate_gmx_itp`](@ref) [`generate_gmx_files`](@ref)

    # Examples
    ```
    julia> ProtoSyn.Peptides.GMX.generate_gmx_topology(pose)
    ✓ All necessary packages were found!
     (...)
    ```
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
        include_atomtypes::Bool    = true)
    
        # Attempts to identify protein and non-protein molecules based on
        # ProteinSelection
        # Generate .itp for protein molecule
        N_protein_molecules = count(ProteinSelection()(pose))
        if N_protein_molecules > 0
            generate_gmx_itp(pose, ProteinSelection(),
                molecule_itp_filename  = protein_itp_filename,
                overwrite          = overwrite,
                keep_temp_files    = keep_temp_files,
                gmx_histidine_type = gmx_histidine_type,
                gmx_water_model    = gmx_water_model,
                gmx_forcefield     = gmx_forcefield)
        end

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
            Base.write(topol, "; Topology file generated by ProtoSyn.jl\n\n")
            Base.write(topol, "; Include forcefield parameters\n")
            Base.write(topol, "#include \"$(gmx_forcefield).ff/forcefield.itp\"\n")
            # println("N_non_protein_molecules: $(N_non_protein_molecules) ")
            if N_non_protein_molecules > 0 && include_atomtypes
                Base.write(topol, "#include \"$(atomtypes_itp_filename)\"\n")
            end
            Base.write(topol, "\n; Include water topology\n")
            Base.write(topol, "#include \"$(gmx_forcefield).ff/$(gmx_water_model).itp\"\n")
            Base.write(topol, "\n; Include topology for ions\n")
            Base.write(topol, "#include \"$(gmx_forcefield).ff/ions.itp\"\n")
            if N_protein_molecules > 0
                Base.write(topol, "\n; Include topology for protein\n")
                Base.write(topol, "#include \"$(protein_itp_filename)\"\n")
            end
            if N_non_protein_molecules > 0
                Base.write(topol, "\n; Include topology for ligands\n")
                Base.write(topol, "#include \"$(molecule_itp_filename)\"\n")
            end
            Base.write(topol, "\n[ system ]\n; Name\n$(pose.graph.name)\n")
            Base.write(topol, "\n[ molecules ]\n; Compound        #mols\n")
            if N_protein_molecules > 0
                Base.write(topol, @sprintf("%-19s %-d\n", protein_itp_filename[1:(end-4)], 1))
            end
            if N_non_protein_molecules > 0
                Base.write(topol, @sprintf("%-19s %-d\n", molecule_itp_filename[1:(end-4)], 1))
            end
        end

        if !keep_temp_files
            isfile("conf.gro")  && rm("conf.gro")
            isfile("posre.itp") && rm("posre.itp")
        end
    end


    struct GMXSA <: ProtoSyn.SearchAlgorithm end

    function (sa::GMXSA)(atom::Atom, stack::Vector{Atom})

        # In THR, CG2 comes before OG1
        if atom.container.name == "THR" && atom.name == "CB"
            bonds = copy(ProtoSyn.sort_children(atom))

            if length(bonds) === 2 # Case no hydrogens
                bonds = bonds[[2, 1]]
            elseif length(bonds) === 3 # Normal case
                bonds = bonds[[1, 3, 2]]
            end
        elseif atom.container.name == "CYS" && atom.name == "CB"
            bonds = copy(ProtoSyn.sort_children(atom))
            S = [a for a in bonds if a.name == "SG"][1]

            if length([a for a in S.bonds if a.symbol == "H"]) === 0 # Case there's a sulfide bond
                bonds = bonds[[2, 3, 1]]
            end
        else
            bonds = copy(ProtoSyn.sort_children(atom))
        end

        if atom.name == "CA"
            c_index = findfirst((atom) -> atom.name == "C", bonds)
            if c_index === nothing
                @warn "Tried to sort $atom children following IUPAC convention, but no \"C\" atom was found."
            else
                push!(bonds, bonds[c_index])
                deleteat!(bonds, c_index)
            end
        end

        return vcat(bonds, stack)
    end

    """
        (ProtoSyn.Peptides.GMX.GMXS)(atom::Atom, stack::Vector{Atom})

    Gromacs-like algorithm for [`travel_graph`](@ref). Correctly sorts the
    given [`Atom`](@ref) `atom` children instances and concatenates with the current
    `stack` as expected by the Gromacs external package.

    # Examples
    ```
    julia> ProtoSyn.Peptides.GMX.GMXS(pose.graph[1, 1, 1], Vector{Atom}())
    4-element Vector{Atom}:
     Atom{/2a3d:34300/A:1/MET:1/H1:2}
     Atom{/2a3d:34300/A:1/MET:1/H2:3}
     Atom{/2a3d:34300/A:1/MET:1/H3:4}
     Atom{/2a3d:34300/A:1/MET:1/CA:5}
    ```
    """
    GMXS = GMXSA()

    """
        sort_atoms_and_graph_gmx(pose::Pose)

    Sorts [`Atom`](@ref) instances in the encompassing `AbstractContainer`
    structures and the resulting [Graph](@ref graph-types) for Gromacs simulations (Gromacs
    expects a given atom order to match with the `gmx pdb2gmx` and `acpype`
    generated .itp and .top files - this does not necessarilly match the IUPAC
    conventions).

    # See also
    [`assign_gmx_atom_names!`](@ref) [`generate_gmx_files`](@ref)

    # Examples
    ```
    julia> ProtoSyn.Peptides.GMX.sort_atoms_and_graph_gmx(pose)
    Pose{Topology}(Topology{/2a3d:60945}, State{Float64}:
     Size: 1140
     i2c: false | c2i: false
     Energy: Dict(:Total => Inf)
    )
    ```
    """
    function sort_atoms_and_graph_gmx(pose::Pose)

        _pose = copy(pose)

        # Deal with Prolines
        for container in rn"PRO"(_pose, gather = true)

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

        # Deal with Tryptophans
        for container in rn"TRP"(_pose, gather = true)

            @assert container.name == "TRP" "Can't sort Tryptophan atoms and graph on a non-Tryptophan residue ($container)."

            # Manually identify necessary heavy-atoms
            NE1 = container["NE1"]; NE1 === nothing && error("Tried to sort atoms and graph in residue $container, but no \"NE1\" atom was found. Check atom names: consider using Peptides.assign_default_atom_names!.")
            CE2 = container["CE2"]; CE2 === nothing && error("Tried to sort atoms and graph in residue $container, but no \"CE2\" atom was found. Check atom names: consider using Peptides.assign_default_atom_names!.")
            CH2 = container["CH2"]; CH2 === nothing && error("Tried to sort atoms and graph in residue $container, but no \"CH2\" atom was found. Check atom names: consider using Peptides.assign_default_atom_names!.")
            CZ3 = container["CZ3"]; CZ3 === nothing && error("Tried to sort atoms and graph in residue $container, but no \"CZ3\" atom was found. Check atom names: consider using Peptides.assign_default_atom_names!.")
            CE3 = container["CE3"]; CE3 === nothing && error("Tried to sort atoms and graph in residue $container, but no \"CE3\" atom was found. Check atom names: consider using Peptides.assign_default_atom_names!.")
            CD2 = container["CD2"]; CD2 === nothing && error("Tried to sort atoms and graph in residue $container, but no \"CD2\" atom was found. Check atom names: consider using Peptides.assign_default_atom_names!.")

            # Manually assign parents
            ProtoSyn.popparent!(CE2); ProtoSyn.setparent!(CE2, NE1)
            ProtoSyn.popparent!(CZ3); ProtoSyn.setparent!(CZ3, CH2)
            ProtoSyn.popparent!(CE3); ProtoSyn.setparent!(CE3, CZ3)
            ProtoSyn.popparent!(CD2); ProtoSyn.setparent!(CD2, CE3)    
        end
        
        # Deal with other aromatic residues
        aromatic_residues = rn"TYR|PHE|HIS|HIE"r(_pose, gather = true)

        for residue in aromatic_residues
            p = residue["N"].parent
            ProtoSyn.infer_parenthood!(residue, overwrite = true, start = residue["N"], linear_aromatics = true)
            ProtoSyn.popparent!(residue["N"])
            ProtoSyn.setparent!(residue["N"], p)
        end

        # Sort atoms in the container according to the new parenthoods
        # The GMXSA already deals with threonines
        ProtoSyn.sort_atoms_by_graph!(_pose, search_algorithm = GMXSA())

        # Re-calculate internal coordinates based on new parenthoods
        ProtoSyn.request_c2i!(_pose.state)
        sync!(_pose)

        return _pose
    end


    """
        assign_gmx_atom_names!(pose::Pose, selection::Opt{AbstractSelection} = nothing)

    Rename [`Atom`](@ref) instances in the given [`Pose`](@ref) `pose` to match
    with the `gmx pdb2gmx` and `acpype` generated .itp and .top files - this
    does not necessarilly match the IUPAC conventions. If an `AbstractSelection`
    `selection` is provided, consider only the selected [`Atom`](@ref) instances
    (any given `selection` will be promoted to be of [`Atom`](@ref) type - see
    [`ProtoSyn.promote`](@ref)).

    # See also
    [`assign_default_atom_names!`](@ref ProtoSyn.Peptides.assign_default_atom_names!) [`rename!`](@ref ProtoSyn.rename!) [`sort_atoms_and_graph_gmx`](@ref) [`generate_gmx_files`](@ref)
    
    # Examples
    ```
    julia> ProtoSyn.Peptides.GMX.assign_gmx_atom_names!(pose)
    Pose{Topology}(Topology{/2a3d:55318}, State{Float64}:
     Size: 1140
     i2c: false | c2i: false
     Energy: Dict(:Total => Inf)
    )
    ```
    """
    function assign_gmx_atom_names!(pose::Pose, selection::Opt{AbstractSelection} = nothing)

        ProtoSyn.Peptides.assign_default_atom_names!(pose, selection & ProteinSelection())

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

        ct_sele     = ProtoSyn.promote((sele & DownstreamTerminalSelection{Residue}()), Residue)
        c_terminals = ct_sele(pose, gather = true)
        for c_terminal in c_terminals
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


    """
        generate_gmx_files(pose::Pose; protein_itp_filename::String ="protein.itp", output_pdb_filename::String ="system.pdb", atomtypes_itp_filename::String ="atomtypes.itp", molecule_itp_filename::String ="molecule.itp", overwrite::Bool =true, keep_temp_files::Bool =false, gmx_histidine_type::Int =0, gmx_water_model::String ="tip3p", gmx_forcefield::String ="amber99sb-ildn", topology_filename::String ="topol.top", include_atomtypes::Bool =false)
    
    Generate all the necessary files to launch a Gromacs MD simulation (.top,
    .itp and correctly ordered/named .pdb files). For topology generation,
    please check [`generate_gmx_topology`](@ref) documentation: all arguments
    used in that method are also input arguments of this one. For .pdb
    generation, ProtoSyn attempts to sort (using
    [`sort_atoms_and_graph_gmx`](@ref)) and rename (using
    [`assign_gmx_atom_names!`](@ref)) [`Atom`](@ref) instances in accordance to
    the expected values in Gromacs before exporting the given [`Pose`](@ref)
    `pose` to a `output_pdb_filename` file. Check these methods documentation
    for more details.

    # Examples
    ```
    julia> ProtoSyn.Peptides.GMX.generate_gmx_files(pose)
    ✓ All necessary packages were found!
     (...)
    ```
    """
    function generate_gmx_files(pose::Pose;
        selection::Opt{AbstractSelection} = nothing,
        protein_itp_filename::String      = "protein.itp",
        output_pdb_filename::String       = "system.pdb",
        atomtypes_itp_filename::String    = "atomtypes.itp",
        molecule_itp_filename::String     = "molecule.itp",
        overwrite::Bool                   = true,
        keep_temp_files::Bool             = false,
        gmx_histidine_type::Int           = 0,
        gmx_water_model::String           = "tip3p",
        gmx_forcefield::String            = "amber99sb-ildn",
        topology_filename::String         = "topol.top",
        include_atomtypes::Bool           = true)

        p = copy(pose)
        if selection !== nothing
            ProtoSyn.pop_atoms!(p, !selection)
        end
        ProtoSyn.Peptides.cap!(p, ProteinSelection()) # GMX templates always have caps
        ProtoSyn.Peptides.GMX.assign_gmx_atom_names!(p, ProteinSelection()) # Assign names first, as graph changes after
        p = ProtoSyn.Peptides.GMX.sort_atoms_and_graph_gmx(p)
        # print(p.graph[1, 1].items)
        generate_gmx_topology(p,
            protein_itp_filename   = protein_itp_filename,
            atomtypes_itp_filename = atomtypes_itp_filename,
            molecule_itp_filename  = molecule_itp_filename,
            overwrite              = overwrite,
            keep_temp_files        = keep_temp_files,
            gmx_histidine_type     = gmx_histidine_type,
            gmx_water_model        = gmx_water_model,
            gmx_forcefield         = gmx_forcefield,
            topology_filename      = topology_filename,
            include_atomtypes      = include_atomtypes)

        ProtoSyn.write(p, output_pdb_filename)
    end # function


    """
    """
    function configure_mdp(filename::String, steps::Int = 5000, print_every::Int = 500, dt::T = 0.001) where {T <: AbstractFloat}
        content = readlines(filename)
        open(filename, "w") do io
            for line in content
                if startswith(line, "nsteps")
                    line *= "$steps"
                elseif startswith(line, "dt")
                    line *= "$dt"
                elseif startswith(line, "nstxout") | startswith(line, "nstenergy") | startswith(line, "nstlog") | startswith(line, "nstxout-compressed")
                    line *= "$print_every"
                end
                write(io, line*"\n")
            end
        end
    end

    """
    # TODO
    """
    function run_minimization(pose::Pose;
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
        include_atomtypes::Bool    = true,
        mdps_dir::String           = "/home/jpereira/scripts/mdps",
        verbose::Bool              = true,
        copy_log::Bool             = false,
        steps::Int                 = 5000,
        print_every::Int           = 500,
        box_threshold::T           = 1.0) where {T <: AbstractFloat}

        wp = pwd()
        d  = mktempdir()
        cd(d)
        generate_gmx_files(pose,
            protein_itp_filename   = protein_itp_filename,
            output_pdb_filename    = output_pdb_filename,
            atomtypes_itp_filename = atomtypes_itp_filename,
            molecule_itp_filename  = molecule_itp_filename,
            overwrite              = overwrite,
            keep_temp_files        = keep_temp_files,
            gmx_histidine_type     = gmx_histidine_type,
            gmx_water_model        = gmx_water_model,
            gmx_forcefield         = gmx_forcefield,
            topology_filename      = topology_filename,
            include_atomtypes      = include_atomtypes)

        mol_size = maximum(ProtoSyn.Calculators.distance_matrix(pose)) * (T(0.2) * (T(1.0) + box_threshold / (T(2.0)))) # 200% + threshold% the maximum size
        name     = output_pdb_filename[1:(end-4)]
        output_file_pdb_filename_box = name*"_box.pdb"
        ProtoSyn.GMX.add_bounding_box(output_pdb_filename, output_file_pdb_filename_box, mol_size)
        mdp = Base.joinpath(mdps_dir, "minimization_1.mdp")
        cp(mdp, pwd()*"/minimization_1.mdp", force = true)
        configure_mdp("minimization_1.mdp", steps, print_every, 0.001)
        tpr = name*".tpr"
        trr = name*".trr"
        v   = verbose ? nothing : devnull
        redirect_stdio(stderr = v, stdout = v) do
            run(`gmx grompp -f minimization_1.mdp -c $output_file_pdb_filename_box -p $topology_filename -o $tpr -maxwarn 10`)
        end
        redirect_stdio(stderr = v, stdout = v) do
            run(`gmx mdrun -deffnm $name -s $tpr`)
        end

        # Get number of steps
        function count_steps(filename::String)
            open(filename, "r") do io
                for line in eachline(io)
                    if startswith(line, "Steepest Descents converged")
                        elems = split(line)
                        if elems[4] === "machine"
                            return parse(Int, elems[7])
                        else
                            return parse(Int, elems[8])
                        end
                    end
                end

                return nothing
            end
        end

        N_conv = count_steps("$name.log")
        N_conv = N_conv !== nothing ? N_conv : steps

        cp(tpr, wp*"/$tpr", force = true)
        cp(trr, wp*"/$trr", force = true)
        if copy_log
            cp(name*".log", wp*"/$name.log", force = true)
        end
        cd(wp)

        redirect_stdio(stderr = v, stdout = v) do
            run(pipeline(`echo 0 0`, `gmx trjconv -f $trr -o $(name*"_last_frame.pdb") -s $tpr -pbc mol -center -b $N_conv`))
        end
    end

    include("gromacs-launch.jl")

end # module