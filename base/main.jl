using ProtoSyn
using Printf
using LinearAlgebra

# ------------------ System Set-up Stage
include("src/config.jl")                                                  # Load configuration file
state, metadata     = Common.load_from_pdb(input_pdb)                     # Load structure and metadata from PDB file
amber_topology      = Forcefield.Amber.load_from_json(input_amber_json)   # Load Amber topology from JSON file
Common.apply_ss!(state, metadata, ss)                                     # Apply SS from string in the configuration file
contact_restraints  = Forcefield.Restraints.load_distance_restraints_from_file(
    input_contact_map, metadata, k = contact_force_constant,
    threshold = contact_threshold, min_distance = contact_min_distance)   # Load contact restraints from contacts file
dihedral_restraints = Forcefield.Restraints.lock_block_bb(
    metadata, fbw = dihedral_fb_width, k = dihedral_force_constant)       # Lock blocks based on the defined secondary structure
nb_dhs              = filter(x -> x.residue.ss == Common.SS.COIL, metadata.dihedrals)       # Define unlocked dihedrals
nb_phi_dhs          = filter(x -> x.dtype == Common.DIHEDRAL.phi, nb_dhs)                   # Define uncloked PHI dihedrals
solv_pairs          = Forcefield.CoarseGrain.compile_solv_pairs(metadata.dihedrals, λ_eSol) # Compile solvation pairs from dihedrals
hb_groups           = Forcefield.CoarseGrain.compile_hb_groups(metadata.atoms, λ_eH)        # Compile hydrogen bonding groups
include("src/callbacks.jl")                                               # Load callbacks file
include("src/drivers.jl")                                                 # Load drivers file

ref_state, ref_metadata = Common.load_from_pdb(ref_native)                # Load reference structure and metadata from PDB file
amber_cr_cg_evaluator_nc!(ref_state, false)                                  # Calculate reference energy
_print_energy_components(-1, ref_state, initial_minimizer, 0.0, "(TRGT)") # Print reference energy to file

# ------------------ Production Stage
initial_minimizer.run!(state, initial_minimizer)                          # Initial minimization
search_driver.run!(state, search_driver)                                  # ILSRR algorithm main run
printstyled(@sprintf("%s\n", "."^150), color=:grey)