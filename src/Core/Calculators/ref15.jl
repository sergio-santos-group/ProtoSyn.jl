module REF15

    using ProtoSyn
    using ProtoSyn.Calculators: EnergyFunctionComponent
    using PyCall

    const pyrosetta = PyNULL()
    const ref15     = PyNULL()

    function __init__()
        try
            copy!(pyrosetta, pyimport("pyrosetta"))
        catch LoadError
            if ENV["JULIA_PROTOSYN_WARN_NON_AVALIABLE_EFC"] === "true"
                println()
                @warn """
                📍 ProtoSyn was not able to identify `pyrosetta` in this system.
                PyCall is currently configured to use the Python version at $(PyCall.current_python()).
                In order to use the REF15 energy function component, make sure to:
                    - Set ENV["PYTHON"] to the path of python executable you wish to use, run Pkg.build("PyCall") and re-launch Julia and ProtoSyn.
                    - Make sure `pyrosetta` is installed in the machine trying to load ProtoSyn.

                In order to install `pyrosetta`, follow the following instructions:
                (1) Obtain a free academic license (https://els2.comotion.uw.edu/product/pyrosetta)
                (2) Add the pyrosetta channel to conda, by adding it to ~/.condarc. Be sure to change the USERNAME and PASSWORD to the provided credentials in the license.
                    channels:
                    - https://USERNAME:PASSWORD@conda.graylab.jhu.edu
                    - conda-forge
                    - defaults
                (3) Install pyrosetta: conda install pyrosetta
                (4) Re-launch Julia and ProtoSyn

                ProtoSyn will continue loading, but the `Calculators.REF15` module will be unavailable.
                To surpress further warnings for unavailable energy function componenets, set the JULIA_PROTOSYN_WARN_NON_AVALIABLE_EFC environment flag and re-launch Julia and ProtoSyn. 
                \$ export JULIA_PROTOSYN_WARN_NON_AVALIABLE_EFC=false
                Optionally, add the above line to ~/.bashrc to persistently supress warnings in further sessions.

                """
            end
        end

        if pyrosetta != PyNULL()
            # Init PyRosetta -> redirect output to a devnull stdout to prevent
            # clutering of the terminal with copyright information
            sys = pyimport("sys")
            os  = pyimport("os")
            sys.stdout = open(os.devnull, "w")
            pyrosetta.init(extra_options = "-mute all")
            sys.stdout = sys.__stdout__

            # Set the REF15 score function
            # `get_score_function` takes an optional Bool argument, in this case
            # `is_full_atom = True`
            copy!(ref15, pyrosetta.get_score_function(true))
        end
    end


    """
    # TODO DOCUMENTATION
    Uses the file system
    Does not care for acceleration type
    Selection?
    """
    function calc_ref15(::Type{<: ProtoSyn.AbstractAccelerationType}, pose::Pose, selection::Opt{AbstractSelection}, update_forces::Bool = false; rosetta_pose::Opt{PyCall.PyObject})

        # * Another (faster) option is to create a pyrosetta poe without using
        # * the filesystem. For example, using the pose_from_sequence method.
        # * This, however, had to take into consideration the terminal cap
        # * status, presence of tautomers, hydrogentation status, missing atoms
        # * in either ProtoSyn or PyRosetta poses, etc.

        # * For some unknown reason, when using a pre-existing rosetta_pose, the
        # * resulting energy is slightly different.
        
        if rosetta_pose === nothing
            # A.1. Write pose to a temporary file. File is added to the tempdir()
            # directory and deleted after completion of the current Julia process
            filename, io = Base.Filesystem.mktemp()
            ProtoSyn.write(io, pose.graph, pose.state, ProtoSyn.PDB)
            close(io)

            # A.2. Load pose as a pyrosetta pose
            _pose = pyrosetta.pose_from_pdb(filename)
        else
            # B.1. Update atomic coordinates on the pre-existing rosetta_pose
            _pose = rosetta_pose
            i = 0
            for r in _pose.residues
                for (atom_index, atom) in enumerate(r.atoms())
                    r.is_virtual(atom_index) && begin
                        continue
                    end
                    i += 1
                    xyz = pyrosetta.rosetta.numeric.xyzVector_double_t(pose.state[i].t...)
                    r.set_xyz(atom_index, xyz)
                end
            end
        end

        _pose.update_residue_neighbors()
        _pose.dump_pdb("teste.pdb")

        # 3. Return the evaluated energy
        return ref15(_pose), nothing
    end

    calc_ref15(pose::Pose, selection::Opt{AbstractSelection}, update_forces::Bool = false) = begin
        calc_ref15(ProtoSyn.acceleration.active, pose, selection, update_forces)
    end

    """
    """
    function fixate_rosetta_pose!(efc::EnergyFunctionComponent{T}, pose::Pose) where {T <: AbstractFloat}
        if (:rosetta_pose in keys(efc.settings))
            if efc.settings[:rosetta_pose] === PyCall.PyObject
                @warn "It seems the current EnergyFunctionComponent already has a set :rosetta_pose setting. ProtoSyn will replace the pre-exisitng :rosetta_pose."
            end

            # 1. Write pose to a temporary file. File is added to the tempdir()
            # directory and deleted after completion of the current Julia process
            filename, io = Base.Filesystem.mktemp()
            ProtoSyn.write(io, pose.graph, pose.state, ProtoSyn.PDB)
            close(io)

            # 2. Load pose as a pyrosetta pose
            _pose = pyrosetta.pose_from_pdb(filename)

            # 3. Save as the efc setting
            efc.settings[:rosetta_pose] = _pose
        end
    end

    """
    """
    function get_default_ref15(;α::T = 1.0) where {T <: AbstractFloat}
        return EnergyFunctionComponent(
            "REF15",
            calc_ref15,
            nothing,
            Dict{Symbol, Any}(:rosetta_pose => nothing),
            α,
            true)
    end
end