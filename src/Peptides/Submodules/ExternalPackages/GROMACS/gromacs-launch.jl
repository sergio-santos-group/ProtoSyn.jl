mutable struct GMX_Step{T <: AbstractFloat}
    name::String
    time_step::T
    n_steps::Int
    print_every::Int
    overwrite::Bool
end

function Base.show(io::IO, gmxs::GMX_Step{T}, level_code::Opt{LevelCode} = nothing) where {T <: AbstractFloat}
    init_level_code = level_code === nothing ? LevelCode() : level_code
    init_lead       = ProtoSyn.get_lead(level_code)

    println(io, init_lead*"∼  GMX Simulation Step:")

    lead       = ProtoSyn.get_lead(level_code)
    inner_lead = ProtoSyn.get_inner_lead(level_code)

    println(io, inner_lead*"+"*repeat("-", 51)*"+")
    @printf(io, "%s|%-15s | %-30s   |\n", inner_lead, " Name", "$(gmxs.name)")
    println(io, inner_lead*"+"*repeat("-", 51)*"+")
    @printf(io, "%s|%-15s | %-30s   |\n", inner_lead, " Time step", "$(gmxs.time_step) ps")
    @printf(io, "%s|%-15s | %-30s   |\n", inner_lead, " Nº of steps", "$(gmxs.n_steps)")
    time_scale = gmxs.n_steps * gmxs.time_step
    time_scale_units = "ps"
    if time_scale >= 1000.0
        time_scale_units = "ns"
        time_scale = time_scale / 1000.0
    end
    @printf(io, "%s|%-15s | %-30s   |\n", inner_lead, " Time scale", @sprintf("%.2f %s", time_scale, time_scale_units))
    @printf(io, "%s|%-15s | %-30s   |\n", inner_lead, " Print every", "$(gmxs.print_every)")
    @printf(io, "%s|%-15s | %-30s   |\n", inner_lead, " Overwrite", "$(gmxs.overwrite)")
    println(io, inner_lead*"+"*repeat("-", 51)*"+")
end

# Default GMX steps
minimization_1    = GMX_Step(     "1-minimization", 0.002,   50_000,  1_000, true)
minimization_2    = GMX_Step(     "2-minimization", 0.002,   50_000,  1_000, true)
nvt_equilibration = GMX_Step("3-nvt-equilibration", 0.002,  250_000, 10_000, true)
npt_equilibration = GMX_Step("4-npt-equilibration", 0.002,  250_000, 10_000, true)
collection        = GMX_Step(       "5-collection", 0.002, 5000_000, 10_000, true)
single_point      = GMX_Step(      "0-single-step", 0.000,        0,      1, true) # Note: Print every can't be 0.0


"""
# TODO Documentation
"""
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


"""
# TODO Documentation
"""
function launch_md_simulation(pose::Pose;
    gmx_steps::Vector{GMX_Step}       = Vector{GMX_Step}(
        [minimization_1,
        minimization_2,
        nvt_equilibration,
        npt_equilibration,
        collection]),
    selection::Opt{AbstractSelection} = nothing,
    mdps_dir::String                  = "/home/jpereira/scripts/mdps",
    simulation_name::Opt{String}      = nothing,
    overwrite::Bool                   = false,
    protein_itp_filename::String      = "protein.itp",
    gmx_histidine_type::Int           = 0,
    box_threshold::T                  = 1.0,
    include_atomtypes::Bool           = true,
    verbose::Bool                     = true,
    add_solvent::Bool                 = true) where {T <: AbstractFloat}

    function check_file(filename::String)
        file_exists = false
        if isfile(filename)
            printstyled("File $filename was found in the current working directory.\n", color = :cyan)
            if overwrite
                printstyled("Overwrite flag is set to $overwrite - Removing pre-existing simulation data.\n", color = :yellow)
            else
                printstyled("Overwrite flag is set to $overwrite - Maintaining pre-existing simulation data.\n", color = :green)
                file_exists = true
            end
        end

        return file_exists
    end

    # 0. Define defaults
    cwd = pwd()
    if simulation_name === nothing
        simulation_name = pose.graph.name*"_md"
    end
    filename_2 = simulation_name*"_pre_sol.pdb"
    filename_3 = simulation_name*"_pre_ion.pdb"
    filename_4 = simulation_name*"_pre_min.pdb"

    # 1. Check the existence of pre-existing simulation folders
    if isdir(simulation_name)
        printstyled(" | Simulation folder $simulation_name was found in the current working directory.\n", color = :cyan)
        if overwrite
            printstyled("   Overwrite flag is set to $overwrite - Removing pre-existing simulation data.\n", color = :yellow)
            rm(simulation_name, force = true, recursive = true)
        else
            printstyled("   Overwrite flag is set to $overwrite - Maintaining pre-existing simulation data.\n", color = :green)
        end
    end

    mkpath(simulation_name)
    cd(simulation_name)
    @info "Now working on $simulation_name"

    # 2. Generate .itp files for system topology
    protein_itp_filename_exists = check_file(protein_itp_filename)

    if !protein_itp_filename_exists
        generate_gmx_files(pose;
            selection = selection,
            protein_itp_filename = protein_itp_filename,
            gmx_histidine_type   = gmx_histidine_type,
            include_atomtypes    = include_atomtypes)
    end

    # 3. Add bounding box
    filename_2_exists = check_file(filename_2)

    if !filename_2_exists
        mol_size = maximum(ProtoSyn.Calculators.distance_matrix(pose)) * (T(0.2) * (T(1.0) + box_threshold / (T(2.0)))) # 200% + threshold% the maximum size
        ProtoSyn.GMX.add_bounding_box("system.pdb", filename_2, mol_size)
    end
    
    # 4. Add solvent
    if add_solvent
        filename_3_exists = check_file(filename_3)

        if !filename_3_exists
            ProtoSyn.GMX.add_solvent(filename_2, filename_3)
        end

        # 5. Add ions
        filename_4_exists = check_file(filename_4)

        if !filename_4_exists
            ion_substitution_group = 13
            if count((!ProteinSelection())(pose)) > 0
                ion_substitution_group = 15
            end
            ProtoSyn.GMX.add_ions(filename_3, filename_4, attempt_auto = ion_substitution_group)
        end
    else
        cp(filename_2, filename_4)
    end

    # 6. Generate restraints
    ProtoSyn.GMX.generate_restraints(filename_4, attempt_auto = 1)
    previous_filename = "../$filename_4"

    # 7. Consume GMX steps
    for step in gmx_steps

        # 7.1 Verify the existence of a simulation folder with the step's name
        data_exists = false
        if isdir(step.name)
            printstyled(" | Simulation folder $(step.name) was found in the current working directory.\n", color = :cyan)
            if step.overwrite
                printstyled("   Overwrite flag is set to $(step.overwrite) - Removing pre-existing simulation data.\n", color = :yellow)
                rm(step.name, force = true, recursive = true)
            else
                printstyled("   Overwrite flag is set to $(step.overwrite) - Maintaining pre-existing simulation data.\n", color = :green)
                data_exists       = true
                previous_filename = "../$(step.name)/$(step.name)_last_frame.pdb"
            end
        end

        if data_exists
            continue
        end

        mkpath(step.name)
        cd(step.name)

        # 7.2 Set-up .mdp file
        mdp = Base.joinpath(mdps_dir, "$(step.name).mdp")
        cp(mdp, Base.joinpath(pwd(), "$(step.name).mdp"), force = true)
        configure_mdp("$(step.name).mdp", step.n_steps, step.print_every, step.time_step)

        # 7.3 Launch simulation step
        tpr = step.name*".tpr"
        trr = step.name*".trr"
        v   = verbose ? nothing : devnull
        redirect_stdio(stderr = v, stdout = v) do
            run(`gmx grompp -f $(step.name).mdp -c $previous_filename -p ../topol.top -o $tpr -maxwarn 10 -r $previous_filename`)
        end
        redirect_stdio(stderr = v, stdout = v) do
            run(`gmx mdrun -deffnm $(step.name) -s $tpr`)
        end

        # 7.4 Extract last PDB frame
        N_conv = count_steps("$(step.name).log") # Only works in minimizations
        N_conv = N_conv !== nothing ? N_conv : convert(Int, step.n_steps / step.print_every)

        redirect_stdio(stderr = v, stdout = v) do
            run(pipeline(`echo 0 0`, `gmx trjconv -f $trr -o $(step.name)_last_frame.pdb -s $tpr -pbc mol -center -b $N_conv -conect`))
        end

        previous_filename = "../$(step.name)/$(step.name)_last_frame.pdb"
        cd("..")
    end

    cd(cwd)
end