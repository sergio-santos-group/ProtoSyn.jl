using ProtoSyn
using Printf
using LinearAlgebra

include("ilsrr_config.jl")
include("ilsrr_tests.jl")

state, metadata        = Common.load_from_pdb(input_pdb)
amber_topology         = Forcefield.Amber.load_from_json(input_amber_json)
Common.fix_proline!(state, metadata.dihedrals)
Common.apply_ss!(state, metadata, ss)
contact_restraints     = Forcefield.Restraints.load_distance_restraints_from_file(input_contact_map, metadata, k = contact_force_constant, threshold = contact_threshold)
dihedral_restraints    = Forcefield.Restraints.lock_block_bb(metadata, fbw = dihedral_fb_width, k = dihedral_force_constant)

nb_dihedrals           = filter(x -> x.residue.ss == Common.SS.COIL, metadata.dihedrals)
nb_phi_dihedrals       = filter(x -> x.dtype == Common.DIHEDRAL.phi, nb_dihedrals) # For crankshaft movements
dihedral_mutator       = Mutators.Dihedral.DihedralMutator(nb_dihedrals, () -> (randn() * dihedral_mutator.step_size), dihedral_p_mut, dihedral_step_size)
crankshaft_mutator     = Mutators.Crankshaft.CrankshaftMutator(nb_phi_dihedrals, () -> (randn() * crankshaft_mutator.step_size), crankshaft_p_mut, crankshaft_step_size)
dihedral_perturbator   = Mutators.Dihedral.DihedralMutator(nb_dihedrals, () -> (randn() * dihedral_mutator.step_size), dihedral_p_mut_prt, dihedral_step_size_prt)
crankshaft_perturbator = Mutators.Crankshaft.CrankshaftMutator(nb_phi_dihedrals, () -> (randn() * crankshaft_mutator.step_size), crankshaft_p_mut_prt, crankshaft_step_size_prt)

include("ilsrr_callbacks.jl")

function my_evaluator!(st::Common.State, do_forces::Bool)
    fill!(st.forces, zero(Float64))
    amber_energy  = Forcefield.Amber.evaluate!(amber_topology, st, cut_off=1.2, do_forces=do_forces)
    other_energy  = Forcefield.Restraints.evaluate!(contact_restraints, st, do_forces=do_forces)
    other_energy += Forcefield.Restraints.evaluate!(dihedral_restraints, st, do_forces=do_forces)
    st.energy.comp["other"] = other_energy
    st.energy.eTotal = amber_energy + other_energy
    return amber_energy + other_energy
end

sampler_sd_driver = Drivers.SteepestDescent.Driver(my_evaluator!, min_n_steps, f_tol, max_step, print_status_sd)
initial_sd_driver = Drivers.SteepestDescent.Driver(my_evaluator!, init_min_n_steps, f_tol, max_step, print_status_sd)
function my_sampler!(st::Common.State)
    Mutators.Dihedral.run!(st, dihedral_mutator)
    Mutators.Crankshaft.run!(st, crankshaft_mutator)
    sampler_sd_driver.run!(st, sampler_sd_driver)
end

function my_pertubator!(st::Common.State)
    dm = Mutators.Dihedral.run!(st, dihedral_perturbator)
    cm = Mutators.Crankshaft.run!(st, crankshaft_perturbator)
    println("(ILSRR) Performing structure perturbation ▶️ Dihedral: $dm | Crankshaft: $cm")
    sampler_sd_driver.run!(st, sampler_sd_driver)
end

inner_cycle_driver = Drivers.MonteCarlo.Driver(my_sampler!, my_evaluator!, temperature, n_inner_steps, print_status_mc, quench_temperature, print_structure_ic, adjust_step_size)
ilsrr_driver = Drivers.ILSRR.Driver(inner_cycle_driver, my_evaluator!, my_pertubator!, temperature, n_outer_steps, reset_temperature, reset_step_size, print_outer_best)
# identify_angles_under_stress(metadata.dihedrals, dihedral_restraints)
# Print.as_pdb(xyz_destination, state, metadata)
# initial_sd_driver.run!(state, initial_sd_driver, print_structure_min)
ilsrr_driver.run!(state, ilsrr_driver)
