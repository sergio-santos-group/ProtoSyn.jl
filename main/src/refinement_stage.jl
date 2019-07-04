#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#   O U T E R   D R I V E R :    S I M U L A T E D   A N N E A L I N G
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# ------------------------------- Minimizers ----------------------------------
rfnm_sd_evaluator = Forcefield.Evaluator(components = [amber_top, dihedral_restraints])

rfnm_sd_minimizer = Drivers.SteepestDescent.DriverConfig(
    evaluator   = rfnm_sd_evaluator,
    n_steps     = 2000,
    nblist_freq = 50,
    f_tol       = 0.2,
    max_step    = 0.01)


# ------------------------------- Evaluator -----------------------------------
rfnm_cycle_evaluator = Forcefield.Evaluator(
    components = [amber_top, hb_network_R, contact_restraints_R, solv_pairs_R])


# -------------------------------- Mutators -----------------------------------
rfnm_dihedral_mutator   = Mutators.Dihedral.MutatorConfig(
    dihedrals     = nb_dhs,
    angle_sampler = () -> (randn() * rfnm_dihedral_mutator.step_size),
    p_mut         = 0.025,
    step_size     = π/32)

rfnm_crankshaft_mutator = Mutators.Crankshaft.MutatorConfig(
    dihedrals     = nb_phi_dhs,
    angle_sampler = () -> (randn() * rfnm_crankshaft_mutator.step_size),
    p_mut         = 0.0135,
    step_size     = π/32)

rfnm_blockrot_mutator = Mutators.Blockrot.MutatorConfig(
    blocks                = metadata.blocks,
    angle_sampler         = () -> (randn() * rfnm_blockrot_mutator.rotation_step_size),
    p_mut                 = 0.25,
    rotation_step_size    = π/16, # rad
    translation_step_size = 0.005, # nm
    rot_axis              = :random,
    trans_axis            = :random,
    n_tries               = 100,
    loop_closer           = rfnm_sd_minimizer)


# -------------------------------- Sampler ------------------------------------
rfnm_cycle_sampler = Mutators.Sampler(
    mutators = [rfnm_dihedral_mutator, rfnm_crankshaft_mutator, rfnm_blockrot_mutator])


# ------------------------------- Callbacks -----------------------------------
rfnm_print_status = Common.@callback 1 function (state, dr_state)
    printstyled(print_energy_components("Refinement", dr_state.step, state.energy, dr_state.temperature)*"\n",
        color = :blue)
end

rfnm_log_energy = Common.@callback 10 function (state, dr_state)
    write(log_energy, print_energy_components("      RFNM", dr_state.step, state.energy, dr_state.temperature)*"\n")
    flush(log_energy)
end

# -------------------------------- Driver -------------------------------------
rfnm_n_steps = 100
rfnm_i_temp  = 200.0

function rfnm_adjust_temperature(step::Int64)
    return -(rfnm_i_temp/rfnm_n_steps) * step + rfnm_i_temp
end

rfnm_sa_driver = Drivers.MonteCarlo.DriverConfig(
    sampler = rfnm_cycle_sampler,
    evaluator = rfnm_cycle_evaluator,
    temperature = rfnm_adjust_temperature,
    n_steps = rfnm_n_steps,
    callbacks = [rfnm_print_status, rfnm_log_energy])