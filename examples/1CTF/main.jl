using ProtoSyn
using Printf
using LinearAlgebra

include("src/config.jl")

state, metadata     = Common.load_from_pdb(input_pdb)
amber_topology      = Forcefield.Amber.load_from_json(input_amber_json)
Common.fix_proline!(state, metadata.dihedrals)
Common.apply_ss!(state, metadata, ss)

contact_restraints  = Forcefield.Restraints.load_distance_restraints_from_file(input_contact_map, metadata, k = contact_force_constant, threshold = contact_threshold)
dihedral_restraints = Forcefield.Restraints.lock_block_bb(metadata, fbw = dihedral_fb_width, k = dihedral_force_constant)
nb_dihedrals        = filter(x -> x.residue.ss == Common.SS.COIL, metadata.dihedrals)
nb_phi_dihedrals    = filter(x -> x.dtype == Common.DIHEDRAL.phi, nb_dihedrals)

include("src/callbacks.jl")
include("src/evaluators.jl")
include("src/drivers.jl")

# MAIN
initial_minimizer.run!(state, initial_minimizer, print_structure_min)
for stage in 1:n_search_refine_cycles
    println("STAGE: $stage")
    search_driver.run!(state, search_driver)
    refine_driver.run!(state, refine_driver)
end