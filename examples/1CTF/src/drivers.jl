# EVALUATORS ----
function amber_evaluator!(st::Common.State, do_forces::Bool)::Float64
    fill!(st.forces, zero(Float64))
    return Forcefield.Amber.evaluate!(amber_topology, st, cut_off=1.2, do_forces=do_forces)
end

function amber_cr_evaluator!(st::Common.State, do_forces::Bool)::Float64
    fill!(st.forces, zero(Float64))
    amber_e   = Forcefield.Amber.evaluate!(amber_topology, st, cut_off=1.2, do_forces=do_forces)
    contact_e = Forcefield.Restraints.evaluator!(dihedral_restraints, st, do_forces=do_forces)
    energy    = amber_e + contact_e
    st.energy.eTotal = energy
    return energy
end

function amber_cr_dr_evaluator!(st::Common.State, do_forces::Bool)::Float64
    fill!(st.forces, zero(Float64))
    amber_e    = Forcefield.Amber.evaluate!(amber_topology, st, cut_off=1.2, do_forces=do_forces)
    contact_e  = Forcefield.Restraints.evaluator!(contact_restraints, st, do_forces=do_forces)
    dihedral_e = Forcefield.Restraints.evaluator!(dihedral_restraints, st, do_forces=do_forces)
    energy = amber_e + contact_e + dihedral_e
    st.energy.eTotal = energy
    return energy
end


# MINIMIZERS ----
initial_minimizer = Drivers.SteepestDescent.Driver(amber_evaluator!, n_init_min_steps, f_tol, max_step, print_status)
loop_closer       = Drivers.SteepestDescent.Driver(amber_evaluator!, n_loop_closer_steps, f_tol, max_step, print_status)
sampler_minimizer = Drivers.SteepestDescent.Driver(amber_cr_dr_evaluator!, n_refine_min_steps, f_tol, max_step, print_status)


# MUTATORS ----
# SOFT MUTATORS (For MonteCarlo refinement)
reg_dh_mutator  = Mutators.Dihedral.DihedralMutator(nb_dhs, () -> (randn() * soft_dh_mutator.step_size), soft_dh_p_mut, soft_dh_step_size)
reg_cs_mutator  = Mutators.Crankshaft.CrankshaftMutator(nb_phi_dhs, () -> (randn() * soft_cs_mutator.step_size), soft_cs_p_mut, soft_cs_step_size)
soft_br_mutator = Mutators.Blockrot.BlockrotMutator(metadata.blocks, () -> (randn() * soft_br_mutator.step_size), soft_br_p_mut, soft_br_step_size, soft_br_tries, loop_closer)

# REGULAR MUTATORS (For ILSRR inner cycle search)
reg_dh_mutator  = Mutators.Dihedral.DihedralMutator(nb_dhs, () -> (randn() * reg_dh_mutator.step_size), reg_dh_p_mut, reg_dh_step_size)
reg_cs_mutator  = Mutators.Crankshaft.CrankshaftMutator(nb_phi_dhs, () -> (randn() * reg_cs_mutator.step_size), reg_cs_p_mut, reg_cs_step_size)

# HARD MUTATORS (For ILSRR outer cycle perturbation)
hard_dh_mutator = Mutators.Dihedral.DihedralMutator(nb_dhs, () -> (randn() * hard_dh_mutator.step_size), hard_dh_p_mut, hard_dh_step_size)
hard_cs_mutator = Mutators.Crankshaft.CrankshaftMutator(nb_phi_dhs, () -> (randn() * hard_cs_mutator.step_size), hard_cs_p_mut, hard_cs_step_size)


# SAMPLERS ----
# SOFT SAMPLER (For MonteCarlo refinement)
function dh_cs_br_soft_sampler!(st::Common.State)
    dm = Mutators.Dihedral.run!(st, soft_dihedral_mutator)
    cm = Mutators.Crankshaft.run!(st, soft_crankshaft_mutator)
    bm = Mutators.Blockrot.run!(st, soft_blockrot_mutator)
    sampler_minimizer.run!(st, sampler_minimizer)
    return Dict("d" => dm, "c" => cm, "b" => bm)
end

# REGULAR SAMPLER (For ILSRR inner cycle search)
function dh_cs_reg_sampler!(st::Common.State)
    dm = Mutators.Dihedral.run!(st, regular_dihedral_mutator)
    cm = Mutators.Crankshaft.run!(st, regular_crankshaft_mutator)
    return Dict("d" => dm, "c" => cm)
end

# HARD SAMPLER (For ILSRR outer cycle perturbation)
function dh_cs_hard_sampler!(st::Common.State)
    dm = Mutators.Dihedral.run!(st, hard_dihedral_mutator)
    cm = Mutators.Crankshaft.run!(st, hard_crankshaft_mutator)
    return Dict("d" => dm, "c" => cm)
end


# DRIVERS ----
# SEARCH
i_search_driver = Drivers.MonteCarlo.Driver(dh_cs_reg_sampler!, amber_cr_evaluator!, i_search_temp_init, n_i_search_steps,
    print_status_mc, quench_temperature, print_structure_ic, adjust_step_size)
search_driver   = Drivers.ILSRR.Driver(i_search_driver, amber_cr_evaluator!, dh_cs_hard_sampler!, o_search_temp_static, n_o_search_steps,
    reset_temperature, reset_step_size, print_outer_best)

# REFINEMENT
refine_driver   = Drivers.MonteCarlo.Driver(dh_cs_br_soft_sampler!, amber_cr_evaluator!, temperature, n_inner_steps,
    print_status_mc, print_structure_ic, adjust_step_size)