using ProtoSyn
using Printf
using LinearAlgebra

include("src/config.jl")

state, metadata     = Common.load_from_pdb(input_pdb)
amber_topology      = Forcefield.Amber.load_from_json(input_amber_json)
Common.fix_proline!(state, metadata.dihedrals)
Common.apply_ss!(state, metadata, ss)
# metadata.ss     = Common.compile_ss(metadata.dihedrals, ss)
# metadata.blocks = Common.compile_blocks(metadata.residues, ss)

contact_restraints  = Forcefield.Restraints.load_distance_restraints_from_file(input_contact_map, metadata, k = contact_force_constant, threshold = contact_threshold, min_distance = contact_min_distance)
dihedral_restraints = Forcefield.Restraints.lock_block_bb(metadata, fbw = dihedral_fb_width, k = dihedral_force_constant)
nb_dhs              = filter(x -> x.residue.ss == Common.SS.COIL, metadata.dihedrals)
nb_phi_dhs          = filter(x -> x.dtype == Common.DIHEDRAL.phi, nb_dhs)

include("src/callbacks.jl")
include("src/drivers.jl")

# MAIN
Print.as_pdb(outer_best_destination, state, metadata, step = 0)
initial_minimizer.run!(state, initial_minimizer, print_structure_fl)
for cycle in 1:n_search_refine_cycles
    printstyled(@sprintf("\n(%5s) %12s\n", "MAIN", @sprintf("⟳ Cycle: %4d", cycle)), color=:blue)
    printstyled(@sprintf("(%5s) %12s \n%s", "MAIN", "Stage: Random Search", "-"^150), color=:blue)
    search_driver.run!(state, search_driver)
    printstyled(@sprintf("\n(%5s) %12s \n%s\n", "MAIN", "Stage: Refinement", "-"^150), color=:blue)
    refine_driver.run!(state, refine_driver)
end
