using ProtoSyn
using Printf
using LinearAlgebra

include("src/config.jl")

function restart_state!(state::Common.State, input_pdb::String)
    state, metadata = Common.load_from_pdb(input_pdb)
    Common.fix_proline!(state, metadata.dihedrals)
    Common.apply_ss!(state, metadata, ss)
end

state, metadata     = Common.load_from_pdb(input_pdb)
amber_topology      = Forcefield.Amber.load_from_json(input_amber_json)
Common.fix_proline!(state, metadata.dihedrals)
Common.apply_ss!(state, metadata, ss)
# metadata.ss     = Common.compile_ss(metadata.dihedrals, ss)
# metadata.blocks = Common.compile_blocks(metadata.residues, ss)

contact_restraints  = Forcefield.Restraints.load_distance_restraints_from_file(input_contact_map, metadata,
k = contact_force_constant, threshold = contact_threshold, min_distance = contact_min_distance)
dihedral_restraints = Forcefield.Restraints.lock_block_bb(metadata, fbw = dihedral_fb_width, k = dihedral_force_constant)
nb_dhs              = filter(x -> x.residue.ss == Common.SS.COIL, metadata.dihedrals)
nb_phi_dhs          = filter(x -> x.dtype == Common.DIHEDRAL.phi, nb_dhs)

include("src/esol.jl")
include("src/callbacks.jl")
include("src/drivers.jl")

# Define energy target
ref_state, ref_metadata = Common.load_from_pdb(ref_native)
amber_cr_es_evaluator_nc!(ref_state, false)
_print_energy_components(-1, ref_state, initial_minimizer, 0.0, "(TRGT)")

# MAIN
Print.as_pdb(xyz_destination, state, metadata)
initial_minimizer.run!(state, initial_minimizer)
Print.as_pdb(xyz_destination, state, metadata)
for cycle in 1:n_search_refine_cycles
    printstyled(@sprintf("\n(%5s) %12s\n", "MAIN", @sprintf("⟳ Cycle: %4d", cycle)), color=:blue)

    # Search
    printstyled(@sprintf("(%5s) %12s \n%s", "MAIN", "Stage: Random Search", "-"^150), color=:blue)
    search_driver.run!(state, search_driver)

    # Restart
    if cycle < n_search_refine_cycles
        printstyled(@sprintf("\n(%5s) %12s \n%s", "MAIN", "Stage: Restarting ...", "-"^150), color=:blue)
        restart_state!(state, input_pdb)
        initial_minimizer.run!(state, initial_minimizer)
        Print.as_pdb(xyz_destination, state, metadata)
    end
end

printstyled(@sprintf("%s\n", "."^150), color=:grey)