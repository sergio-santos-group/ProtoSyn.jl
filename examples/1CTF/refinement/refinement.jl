using ProtoSyn
using Printf
using LinearAlgebra

include("config.jl")

state, metadata     = Common.load_from_pdb(input_pdb)
amber_topology      = Forcefield.Amber.load_from_json(input_amber_json)
metadata.ss         = Common.compile_ss(metadata.dihedrals, ss)
metadata.blocks     = Common.compile_blocks(metadata.residues, ss)
bb_dhs              = filter(x -> x.dtype <= Common.DIHEDRAL.omega, metadata.dihedrals)

contact_restraints  = Forcefield.Restraints.load_distance_restraints_from_file(input_contact_map, metadata,
    k = contact_force_constant, threshold = contact_threshold, min_distance = contact_min_distance)
dihedral_restraints = Forcefield.Restraints.lock_block_bb(metadata, fbw = dihedral_fb_width, k = dihedral_force_constant)
nb_dhs              = filter(x -> x.residue.ss == Common.SS.COIL, metadata.dihedrals)
nb_phi_dhs          = filter(x -> x.dtype == Common.DIHEDRAL.phi, nb_dhs)
solv_pairs          = Forcefield.CoarseGrain.load_solv_pairs_default(nb_phi_dhs, λ_eSol)

include("../src/callbacks.jl")
include("../src/drivers.jl")

# Aplly sampled backbone to all-atom model
Common.apply_backbone_from_file!(state, bb_dhs, conf_pdb)
initial_minimizer.run!(state, initial_minimizer)
Print.as_pdb(xyz_destination, state, metadata)
amber_cr_dr_es_evaluator_rfnm!(state, false)

# Define energy target
ref_state, ref_metadata = Common.load_from_pdb(ref_native)
amber_cr_dr_es_evaluator_rfnm!(ref_state, false)
_print_energy_components(-1, ref_state, initial_minimizer, 0.0, "(TRGT)")

# Main algorithm
for cycle in 1:n_search_refine_cycles
    printstyled(@sprintf("\n(%5s) %12s\n", "MAIN", @sprintf("⟳ Cycle: %4d", cycle)), color=:blue)
    printstyled(@sprintf("\n(%5s) %12s \n%s\n", "MAIN", "Stage: Refinement", "-"^150), color=:blue)
    refine_driver.run!(state, refine_driver)
end