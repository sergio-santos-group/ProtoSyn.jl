#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                 O U T E R   D R I V E R :    I L S R R
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

include("inner_driver.jl")

# ------------------------------- Minimizers ----------------------------------
outer_sd_evaluator = Forcefield.Evaluator(components = [amber_top, dihedral_restraints])

outer_sd_minimizer = Drivers.SteepestDescent.DriverConfig(
    evaluator   = outer_sd_evaluator,
    n_steps     = 2000,
    nblist_freq = 50,
    f_tol       = 0.2,
    max_step    = 0.01)

# -------------------------- Perturbator Mutators -----------------------------
outer_dihedral_perturbator   = Mutators.Dihedral.MutatorConfig(
    dihedrals     = nb_dhs,
    angle_sampler = () -> (randn() * outer_dihedral_perturbator.step_size),
    p_mut         = 0.2,
    step_size     = π/4)

outer_crankshaft_perturbator = Mutators.Crankshaft.MutatorConfig(
    dihedrals     = nb_phi_dhs,
    angle_sampler = () -> (randn() * outer_crankshaft_perturbator.step_size),
    p_mut         = 0.1,
    step_size     = π/4)

outer_blockrot_mutator = Mutators.Blockrot.MutatorConfig(
    blocks                = metadata.blocks,
    angle_sampler         = () -> (randn() * outer_blockrot_mutator.rotation_step_size),
    p_mut                 = 0.5,
    rotation_step_size    = π/18,
    translation_step_size = 0.05,
    rot_axis              = :random,
    trans_axis            = :random,
    n_tries               = 100,
    loop_closer           = outer_sd_minimizer)

# --------------------------- Perturbator Sampler -----------------------------
smpl_ilsrr_perturbator = Mutators.Sampler(
    mutators = [outer_dihedral_perturbator, outer_crankshaft_perturbator]) # ! Took blockrot out


# ------------------------------- Callbacks -----------------------------------
print_status_ilsrr = Common.@callback 1 function (state, dr_state)
    printstyled(print_energy_keys("ILSRR", dr_state.step)*"\n",
        color = :green)
    printstyled(print_energy_components("ILSRR", dr_state.step, dr_state.best_state.energy, dr_state.temperature)*"\n",
        color = :green)
end

print_smpl_best = Common.@callback 1 function (state, dr_state)
    Print.as_pdb(sampling_bests, dr_state.best_state, metadata)
    flush(sampling_bests)
end

print_status_log = Common.@callback 1 function (state, dr_state)
    write(log_energy, "ILSRR $(dr_state.home_state.energy.total) $(dr_state.best_state.energy.total) \n")
    flush(log_energy)
end

# --------------------------------- Driver ------------------------------------
smpl_ilsrr_driver = Drivers.ILSRR.DriverConfig(
    inner_driver_config = smpl_sa_driver,
    perturbator         = smpl_ilsrr_perturbator,
    temperature         = 1_000.0,
    n_steps             = 10,
    stall_limit         = 5,
    callbacks           = [print_status_ilsrr, print_smpl_best, print_status_log])