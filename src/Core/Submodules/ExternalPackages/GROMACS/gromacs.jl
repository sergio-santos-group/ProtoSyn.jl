"""
# TODO Documentation
"""
module GMX

    using ProtoSyn
    using Printf

    """
    # TODO: Documentation
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
            read(`acpype -v`)
        catch
            error("It seems acpype is not installed in your machine. ProtoSyn tried to call acpype with the following command, but no installation was found: `acpype -v`. Check if acpype installation directory is in this system's PATH.")
        end
    end # function

    """
    # TODO: Documentation
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

end # module