# EVALUATORS ----
# amber- = Amber forcefield parameters
# -cr-   = Contact restraints
# -dr-   = Dihedral restraints
# -cg-   = Coarse-Grained contributions (eH, eSol)
# -nc    = No Coulomb contributions

# Used in the Loop Closer minimizer
function amber_evaluator!(st::Common.State, do_forces::Bool)::Float64
    fill!(st.forces, zero(Float64))
    return Forcefield.Amber.evaluate!(amber_topology, st, cut_off = nonbonded_cut_off, do_forces = do_forces)
end

# Used in the Initial minimizer
function amber_evaluator_nc!(st::Common.State, do_forces::Bool)::Float64
    fill!(st.forces, zero(Float64))
    energy                       = Forcefield.Amber.evaluate!(amber_topology.bonds, state, do_forces = do_forces)
    energy                      += Forcefield.Amber.evaluate!(amber_topology.angles, state, do_forces = do_forces)
    energy                      += Forcefield.Amber.evaluate!(amber_topology.atoms, state, do_forces = do_forces, cut_off = nonbonded_cut_off, eCoulomb_λ = 0.0)
    energy                      += Forcefield.Amber.evaluate!(amber_topology.dihedralsCos, state, do_forces = do_forces)
    st.energy.comp["amber"]      = energy
    st.energy.eTotal             = energy
    return energy
end

# Used in the Sampler minimizer
function amber_cr_dr_evaluator!(st::Common.State, do_forces::Bool)::Float64
    fill!(st.forces, zero(Float64))
    amber_e                      = Forcefield.Amber.evaluate!(amber_topology, st, cut_off=1.2, do_forces=do_forces)
    contact_e                    = Forcefield.Restraints.evaluate!(contact_restraints, st, do_forces=do_forces)
    dihedral_e                   = Forcefield.Restraints.evaluate!(dihedral_restraints, st, do_forces=do_forces)
    energy                       = amber_e + contact_e + dihedral_e
    st.energy.comp["other"]      = contact_e + dihedral_e
    st.energy.eTotal             = energy
    return energy
end

# Main evaluator function for Monte Carlo in Refinement stages
function amber_cr_cg_evaluator!(st::Common.State, do_forces::Bool)::Float64
    fill!(st.forces, zero(Float64))
    amber_e                       = Forcefield.Amber.evaluate!(amber_topology, st, cut_off=1.2, do_forces=do_forces)
    contact_e                     = Forcefield.Restraints.evaluate!(contact_restraints, st, do_forces=do_forces)
    sol_e                         = Forcefield.CoarseGrain.evaluate!(solv_pairs, st, do_forces=do_forces)
    hb_e                          = Forcefield.CoarseGrain.evaluate!(hb_groups, st, do_forces=do_forces)
    st.energy.comp["amber"]       = amber_e
    st.energy.comp["coarseGrain"] = contact_e + sol_e + hb_e
    st.energy.eTotal              = amber_e + contact_e + sol_e + hb_e
    return amber_e + contact_e + sol_e + hb_e
end

# Main evaluator function for Monte Carlo in Search stages
function amber_cr_cg_evaluator_nc!(st::Common.State, do_forces::Bool)::Float64
    fill!(st.forces, zero(Float64))
    amber_e                       = Forcefield.Amber.evaluate!(amber_topology.bonds, st, do_forces = do_forces)
    amber_e                      += Forcefield.Amber.evaluate!(amber_topology.angles, st, do_forces = do_forces)
    amber_e                      += Forcefield.Amber.evaluate!(amber_topology.atoms, st, do_forces = do_forces, cut_off = 1.2, eCoulomb_λ = 0.0)
    amber_e                      += Forcefield.Amber.evaluate!(amber_topology.dihedralsCos, st, do_forces = do_forces)
    contact_e                     = Forcefield.Restraints.evaluate!(contact_restraints, st, do_forces=do_forces)
    sol_e                         = Forcefield.CoarseGrain.evaluate!(solv_pairs, st, do_forces=do_forces)
    hb_e                          = Forcefield.CoarseGrain.evaluate!(hb_groups, st, do_forces=do_forces)
    st.energy.comp["amber"]       = amber_e
    st.energy.comp["coarseGrain"] = contact_e + sol_e + hb_e
    st.energy.eTotal              = amber_e + contact_e + sol_e + hb_e
    return amber_e + contact_e + sol_e + hb_e
end


# MINIMIZERS ----
initial_minimizer = Drivers.SteepestDescent.Driver(amber_evaluator_nc!   , n_init_min_steps   , f_tol, max_step, true, print_status_init) # Used to initialize a new search cycle
loop_closer       = Drivers.SteepestDescent.Driver(amber_evaluator!      , n_loop_closer_steps, f_tol, max_step, true, print_status_loop) # Used to close loops in Blockrot movements
sampler_minimizer = Drivers.SteepestDescent.Driver(amber_cr_dr_evaluator!, n_refine_min_steps , f_tol, max_step, true, print_status_smpl) # Used to perform a small minimization step in Monte Carlo cycles


# MUTATORS ----
# -dh- = Dihedral
# -cs- = Crankshaft
# -br- = Blockrot

# SOFT MUTATORS (For Monte Carlo Refinement)
soft_dh_mutator = Mutators.Dihedral.DihedralMutator(nb_dhs, () -> (randn() * soft_dh_mutator.step_size), soft_dh_p_mut, soft_dh_step_s)
soft_cs_mutator = Mutators.Crankshaft.CrankshaftMutator(nb_phi_dhs, () -> (randn() * soft_cs_mutator.step_size), soft_cs_p_mut, soft_cs_step_s)
soft_br_mutator = Mutators.Blockrot.BlockrotMutator(metadata.blocks, () -> (randn() * soft_br_mutator.step_size), soft_br_p_mut, soft_br_step_s, soft_br_tr_step_s, soft_br_tries, loop_closer)

# REGULAR MUTATORS (For ILSRR inner cycle search)
reg_dh_mutator  = Mutators.Dihedral.DihedralMutator(nb_dhs, () -> (randn() * reg_dh_mutator.step_size), reg_dh_p_mut, reg_dh_step_s)
reg_cs_mutator  = Mutators.Crankshaft.CrankshaftMutator(nb_phi_dhs, () -> (randn() * reg_cs_mutator.step_size), reg_cs_p_mut, reg_cs_step_s)

# HARD MUTATORS (For ILSRR outer cycle PERTURBATOR)
hard_dh_mutator = Mutators.Dihedral.DihedralMutator(nb_dhs, () -> (randn() * hard_dh_mutator.step_size), hard_dh_p_mut, hard_dh_step_s)
hard_cs_mutator = Mutators.Crankshaft.CrankshaftMutator(nb_phi_dhs, () -> (randn() * hard_cs_mutator.step_size), hard_cs_p_mut, hard_cs_step_s)


# SAMPLERS ----
# SOFT SAMPLER (For MonteCarlo refinement)
function dh_cs_br_soft_sampler!(st::Common.State)
    dm = Mutators.Dihedral.run!(st, soft_dh_mutator)
    cm = Mutators.Crankshaft.run!(st, soft_cs_mutator)
    bm = Mutators.Blockrot.run!(st, soft_br_mutator)
    sampler_minimizer.run!(st, sampler_minimizer)
    return Dict("d" => dm, "c" => cm, "b" => bm)
end

# REGULAR SAMPLER (For ILSRR inner cycle search)
function dh_cs_reg_sampler!(st::Common.State)
    dm = Mutators.Dihedral.run!(st, reg_dh_mutator)
    cm = Mutators.Crankshaft.run!(st, reg_cs_mutator)
    # sampler_minimizer.run!(st, sampler_minimizer)
    return Dict("d" => dm, "c" => cm)
end

# HARD SAMPLER (For ILSRR outer cycle perturbation)
function dh_cs_hard_sampler!(st::Common.State)
    dm = Mutators.Dihedral.run!(st, hard_dh_mutator)
    cm = Mutators.Crankshaft.run!(st, hard_cs_mutator)
    printstyled("(ILSRR) Perturbator ⟶\n", color = :blue)
    return Dict("d" => dm, "c" => cm)
end


# DRIVERS ----
# Search Stage
i_search_driver = Drivers.MonteCarlo.Driver(dh_cs_reg_sampler!, amber_cr_cg_evaluator_nc!, i_search_temp_init,
    n_i_search_steps, evaluate_slope_every, evaluate_slope_threshold,true, print_status_mc, quench_temperature,
    print_energy_components, adjust_step_size)

search_driver   = Drivers.ILSRR.Driver(i_search_driver, amber_cr_cg_evaluator_nc!, dh_cs_hard_sampler!,
    o_search_temp_static, n_o_search_steps, continue_after_n_attemps, reset_temperature, reset_step_size,
    print_structure_best, print_structure)

# Refinement Stage
refine_driver   = Drivers.MonteCarlo.Driver(dh_cs_br_soft_sampler!, amber_cr_cg_evaluator!, refine_temp_static,
    n_refine_steps, evaluate_slope_every, evaluate_slope_threshold,true, print_status_refnm, print_structure_best,
    print_structure_rfnm)