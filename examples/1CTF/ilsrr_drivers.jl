# ----------------------
# DEFINE THE DRIVERS:
# ----------------------

# 1. Steepest Descent ---------------------------------------------------------------------
# Load necessary tolopogies
topology = Forcefield.Amber.load_from_json(input_amber_json)
# Define the evaluator
function my_sd_evaluator!(st::Common.State, do_forces::Bool)
    energy = Forcefield.Amber.evaluate!(topology, st, cut_off=1.2, do_forces=do_forces)
    st.energy.comp["amber"] = energy
    st.energy.eTotal = energy
    return energy
end
#Define the driver itself
sd_driver = Drivers.SteepestDescent.SteepestDescentDriver(my_sd_evaluator!, n_steps = min_n_steps)


# 2. Monte Carlo -------------------------------------------------------------------------
# Load necessary topologies
contact_pairs      = Forcefield.Other.load_contact_maps_from_file(input_contact_map, metadata)
# Define the evaluator
function my_evaluator!(st::Common.State, do_forces::Bool)
    Drivers.SteepestDescent.run!(state, sd_driver)
    energy  = Forcefield.Amber.evaluate!(topology, st, cut_off=1.2, do_forces=do_forces)
    st.energy.comp["amber"] = energy
    energy += Forcefield.Other.calc_eContact!(st, contact_pairs, k = 100.0)
    st.energy.eTotal = energy
    return energy
end
# Define the sampler mutators
dihedral_mutator   = Mutators.Dihedral.DihedralMutator(nb_dihedrals, () -> (randn() * dihedral_mutator.step_size), dihedral_p_mut, dihedral_step_size)
crankshaft_mutator = Mutators.Crankshaft.CrankshaftMutator(nb_phi_dihedrals, () -> (randn() * crankshaft_mutator.step_size), crankshaft_p_mut, crankshaft_step_size)
# Define the sampler itself
function my_sampler!(st::Common.State)
    dm = Mutators.Dihedral.run!(st, dihedral_mutator)
    cm = Mutators.Crankshaft.run!(st, crankshaft_mutator)

    return Dict("d" => dm, "c" => cm)
end
# Define the driver itself
mc_driver = Drivers.MonteCarlo.MonteCarloDriver(my_sampler!, my_evaluator!, temperature=temperature, n_steps=n_inner_steps)

# 3. Pertubator -------------------------------------------------------------------------------
# Define the pertubator mutators
dihedral_perturbator = Mutators.Dihedral.DihedralMutator(nb_dihedrals, () -> (rand() * 2 - 1 * dihedral_mutator.step_size), dihedral_p_mut_prt, dihedral_step_size_prt)
crankshaft_perturbator = Mutators.Crankshaft.CrankshaftMutator(nb_phi_dihedrals, () -> (rand() * 2 - 1 * crankshaft_mutator.step_size), crankshaft_p_mut_prt, crankshaft_step_size_prt)
# Define the pertubator itself
function my_perturbator!(st::Common.State)
    dm = Mutators.Dihedral.run!(state, dihedral_perturbator)
    cm = Mutators.Crankshaft.run!(state, crankshaft_perturbator)
    
    return Dict("d" => dm, "c" => cm)
end