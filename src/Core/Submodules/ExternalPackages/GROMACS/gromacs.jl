"""
    GMX (Module)

Holds utilitary methods for preparing and launching GROMCAS simulations.
"""
module GMX

    using ProtoSyn
    using Printf

    """
        check_installation()

    Checks the current machine for necessary packages: `gmx`, `antechamber` &
    `acpype`. Make sure these programs are in the system's PATH.

    # Examples
    ```
    julia> ProtoSyn.GMX.check_installation()
    ✓ All necessary packages were found!
    ```
    """
    function check_installation()

        # Check gmx instalation
        try
            read(`gmx --version`)
        catch
            error("It seems GROMACS is not installed in your machine. ProtoSyn tried to call GROMACS with the following command, but no installation was found: `gmx`. Check if GROMACS installation directory is in this system's PATH.")
        end

        # Check antechamber instalation
        try
            read(`antechamber -L`)
        catch
            error("It seems Antechamber is not installed in your machine. ProtoSyn tried to call Antechamber with the following command, but no installation was found: `antechamber`. Check if Antechamber installation directory is in this system's PATH.")
        end

        # Check acpype instalation
        try
            read(`acpype -h`)
        catch
            error("It seems acpype is not installed in your machine. ProtoSyn tried to call acpype with the following command, but no installation was found: `acpype -v`. Check if acpype installation directory is in this system's PATH.")
        end

        println("✓ All necessary packages were found!")
    end # function


    """
        generate_gmx_itp(pose::Pose, selection::Opt{AbstractSelection}; atomtypes_itp_filename::String = "atomtypes.itp", molecule_itp_filename::String = "molecule.itp", overwrite::Bool = true, keep_temp_files::Bool = false)

    Generate an `.itp` file for the [`Fragment`](@ref) given by applying the
    `AbstractSelection` `selection` on the given [`Pose`](@ref) `pose`. If no
    `selection` is provided, the whole [`Pose`](@ref) `pose` will be considered
    as a single entity for `.itp` generation. This function assumes the selected
    [`Fragment`](@ref) is an unusual organic chemical compound and uses `acpype`
    external software to identify atomtypes and write an `.itp` with AMBER
    forcefield parameters for usage in GROMACS simulations. By default, `acpype`
    prints 2 files: the requested `.itp` (to a file named
    `molecule_itp_filename`) and an `.itp` file with all atomtypes retrieved (to
    a file named `atomtypes_itp_filename`). These atomtypes can sometimes be
    useful to verify the `acpype` automatic attribution or to include in the
    final `topology`. In any of these files are found in the current working
    directory, the function raises an error. If `overwrite` is set to `true`
    (is, by default), the existing files are overwritten instead. Any temporary
    files (starting with "jl_") and directories and deleted after completion,
    unless `keep_temp_files` is set to `true` (`false`, by default). It is
    reccomended that this function is not used to specific organic polymers,
    such as proteins. Specific ProtoSyn modules (such as `Peptides`, in this
    example) should provide more accurate `.itp` generation methods.

    !!! ukw "Note:"
        Other modules, such as Peptides, may provide [`generate_gmx_itp`](@ref) methods for specific molecule types (proteins and peptides, in this example).

    # Examples
    ```
    julia> ProtoSyn.GMX.generate_gmx_itp(pose, rn"CBZ")
    (...)
    ```
    """
    function generate_gmx_itp(pose::Pose, selection::Opt{AbstractSelection};
        atomtypes_itp_filename::String = "atomtypes.itp",
        molecule_itp_filename::String  = "molecule.itp",
        overwrite::Bool            = true,
        keep_temp_files::Bool      = false)

        # Check if necessary packages are installed in this machine
        check_installation()

        # Generate Pose to file from just the selected region. No check on
        # validity on the given selection is done.
        frag      = Pose(fragment(pose, selection))
        temp_path = tempname(".")
        temp_file = temp_path*".pdb"
        ProtoSyn.write(frag, temp_file)

        # Run acpype to generate topology files
        run(`acpype -i $temp_file -a amber`)

        # Separate resulting .itp file in two: [ atomtypes ] + [ moleculetype ]
        function write_atomtypes_and_moleculetype_to_file(atomtypes_itp_filename::String, molecule_itp_filename::String)
            if isfile(atomtypes_itp_filename)
                if overwrite
                    @warn "File $atomtypes_itp_filename already exists in the current directory. ProtoSyn will overwrite existing file. Check if this is the desired behaviour. Consider setting the `overwrite` flag to `false` or setting a new `atomtypes_itp_filename`."
                    open(atomtypes_itp_filename, "w") do io; end
                else
                    error("File $atomtypes_itp_filename already exists in the current directory. Consider setting the `overwrite` flag to `true` or setting a new `atomtypes_itp_filename`.")
                end
            end

            if isfile(molecule_itp_filename)
                if overwrite
                    @warn "File $molecule_itp_filename already exists in the current directory. ProtoSyn will overwrite existing file. Check if this is the desired behaviour. Consider setting the `overwrite` flag to `false` or setting a new `molecule_itp_filename`."
                    open(molecule_itp_filename, "w") do io; end
                else
                    error("File $molecule_itp_filename already exists in the current directory. Consider setting the `overwrite` flag to `true` or setting a new `molecule_itp_filename`.")
                end
            end

            open(atomtypes_itp_filename, "a") do io_atomtypes
                open(molecule_itp_filename, "a") do io_molecule
                    open(joinpath(temp_path*".acpype", temp_path[3:end]*"_GMX.itp"), "r") do io
                        in_atomtypes = false
                        to_skip      = 0
                        for line in eachline(io)
                            if to_skip > 0
                                to_skip -= 1
                                continue
                            end
                            if startswith(line, "[ atomtypes ]")
                                in_atomtypes = true
                            elseif startswith(line, "[ moleculetype ]")
                                in_atomtypes = false
                            end

                            if in_atomtypes
                                Base.write(io_atomtypes, line*"\n")
                            else
                                Base.write(io_molecule, line*"\n")
                                if startswith(line, ";name")
                                    Base.write(io_molecule, @sprintf("%-17s %-d\n", molecule_itp_filename[1:(end-4)], 3))
                                    to_skip += 1
                                end
                            end
                        end
                    end
                end
            end
        end

        write_atomtypes_and_moleculetype_to_file(atomtypes_itp_filename, molecule_itp_filename)

        # Remove .acpype folder
        if !keep_temp_files
            rm(temp_path*".acpype", recursive = true)
            for file in readdir(".")
                endswith(file, "#") && rm(file)
                startswith(file, "jl_") && rm(file)
            end
        end
    end # function


    """
        add_bounding_box(input_filename::String, output_filename::String, size::T, [shape::String = "cubic"], [verbose::Bool = true]) where {T <: AbstractFloat}

    Employ `gmx editconf` to add a bouding box to a given input file
    (`input_filename`), outputs to `output_filename`. The box has the provided
    `shape` and `size` (in nm). If `verbose` is set to `false` (`true`, by
    default), hide the GROMACS output.

    # Examples
    ```
    julia> ProtoSyn.GMX.add_bounding_box("md.pdb", "md_box.pdb", 1.0)
    (...)
    ```
    """
    function add_bounding_box(input_filename::String, output_filename::String,
        size::T, shape::String = "cubic", verbose::Bool = true) where {T <: AbstractFloat}

        v = verbose ? nothing : devnull
        redirect_stdio(stderr = v, stdout = v) do
            run(`gmx editconf -f $input_filename -o $output_filename -c -box $size -bt $shape`)
        end

        return nothing
    end # function


    """
        add_solvent(input_filename::String, output_filename::String, [solvent_type::String = "spc216"], [topol_filename::String = "topol.top"], [verbose::Bool = true])

    Employ `gmx solvate` to add a solvent to a given input file
    (`input_filename`), outputs to `output_filename`. The file contents should
    be incorporated in a bounding box (see [`add_bounding_box`](@ref)). The
    solvent used in given by `solvent_type` (GROMCAS searchs the current working
    directory and the default instalation directory for `.gro` files with this
    name containing structural information of the solvent molecule). Add solvent
    information (such as number of added molecules) to the given
    `topol_filename`. If `verbose` is set to `false` (`true`, by default), hide
    the GROMACS output.

    # Examples
    ```
    julia> ProtoSyn.GMX.add_solvent("md_box.pdb", "md_sol.pdb")
    (...)
    ```
    """
    function add_solvent(input_filename::String, output_filename::String,
        solvent_type::String = "spc216", topol_filename::String = "topol.top", verbose::Bool = true)

        v = verbose ? nothing : devnull
        redirect_stdio(stderr = v, stdout = v) do
            run(`gmx solvate -cp $input_filename -cs $solvent_type.gro -o $output_filename -p $topol_filename`)
        end
        
        return nothing
    end # function


    """
        add_ions(input_filename::String, output_filename::String; [mdp::String = joinpath(ProtoSyn.resource_dir, "ExternalPackages/GROMACS/mdps/ions.mdp")], [topol_filename::String = "topol.top"], [tpr_filename::String = "ions.tpr"], [positive_ion::String = "NA"], [negative_ion::String = "CL"], [number_of_positive_ions::Int = 0], [number_of_negative_ions::Int = 0], [neutral::Bool = true], [attempt_auto::Opt{Int} = nothing], [verbose::Bool = true])
    
        Employ `gmx genion` to add ions to a given input file
    (`input_filename`), outputs to `output_filename`. The file contents should
    be incorporated in a solvated bounding box (see [`add_bounding_box`](@ref)
    and [`add_solvent`](@ref)). Adds `number_of_negative_ions` of `negative_ion`
    and `number_of_positive_ions` of `positive_ion`. Optionally, if `neutral` is
    set to `true`, `gmx genion` will add either `positive_ion` or `negative_ion`
    instances until neutral charge is achieved (reccomended). As a middle step,
    `gmx genion` uses the given ions `mdp` file and writes a temporary
    `tpr_filename`. Add ions information (such as number of added molecules) to
    the given `topol_filename`. If `verbose` is set to `false` (`true`, by
    default), hide the GROMACS output. If an `attempt_auto` is given,
    automatically chooses the SOL group for ion replacement (interactive, by
    default). Note that depending on the system, the SOL group number changes.
    For example, for a protein with ligand, the SOL group should be 15.
    
    # Examples
    ```
    julia> ProtoSyn.GMX.add_ions("md_sol.pdb", "md_ready.pdb", attempt_auto = 15)
    (...)
    ```

    """
    function add_ions(input_filename::String, output_filename::String;
        mdp::String                  = joinpath(ProtoSyn.resource_dir, "ExternalPackages/GROMACS/mdps/ions.mdp"),
        topol_filename::String       = "topol.top",
        tpr_filename::String         = "ions.tpr",
        positive_ion::String         = "NA",
        negative_ion::String         = "CL",
        number_of_positive_ions::Int = 0,
        number_of_negative_ions::Int = 0,
        neutral::Bool                = true,
        attempt_auto::Opt{Int}       = nothing,
        verbose::Bool                = true)

        v = verbose ? nothing : devnull
        charge = neutral ? "-neutral" : "-nn $number_of_negative_ions -np $number_of_positive_ions"
        redirect_stdio(stderr = v, stdout = v) do
            run(`gmx grompp -f $mdp -c $input_filename -p $topol_filename -o $tpr_filename`)
            if attempt_auto !== nothing
                run(pipeline(`echo $attempt_auto`, `gmx genion -s $tpr_filename -p $topol_filename -o $output_filename -pname $positive_ion -nname $negative_ion $charge`))
            else
                run(`gmx genion -s $tpr_filename -p $topol_filename -o $output_filename -pname $positive_ion -nname $negative_ion $charge`)
            end
        end

        return nothing
    end # function

end # module