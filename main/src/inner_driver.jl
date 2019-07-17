#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#   I N N E R   D R I V E R :    S I M U L A T E D   A N N E A L I N G
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ------------------------------- Evaluator -----------------------------------
inner_cycle_evaluator = Forcefield.Evaluator(
    components = [amber_top, hb_network, contact_restraints, solv_pairs])


# -------------------------------- Mutators -----------------------------------
inner_dihedral_mutator   = Mutators.Dihedral.MutatorConfig(
    dihedrals = nb_dhs,
    angle_sampler = () -> (randn() * inner_dihedral_mutator.step_size),
    p_mut = 0.05,
    step_size = π/4)

inner_crankshaft_mutator = Mutators.Crankshaft.MutatorConfig(
    dihedrals = nb_phi_dhs,
    angle_sampler = () -> (randn() * inner_crankshaft_mutator.step_size),
    p_mut = 0.027,
    step_size = π/4)

# Define only Beta sheet blocks
β_blocks = filter(x -> x.type == Common.SS.HELIX, metadata.blocks)

inner_blockrot_mutator = Mutators.Blockrot.MutatorConfig(
    blocks = metadata.blocks, 
    movable_blocks = β_blocks,
    angle_sampler = () -> (randn() * inner_blockrot_mutator.rotation_step_size),
    trans_sampler = () -> (randn() * inner_blockrot_mutator.translation_step_size),
    p_mut = 1.0,
    rotation_step_size = π/32,
    translation_step_size = 0.0,
    rot_axis = :longitudinal,
    n_tries = 100)


# -------------------------------- Sampler ------------------------------------
inner_cycle_sampler = Mutators.Sampler(
    mutators = [inner_dihedral_mutator, inner_crankshaft_mutator, inner_blockrot_mutator])


# ------------------------------- Callbacks -----------------------------------
print_status_sa = Common.@callback 5000 function (state, dr_state)
    println(print_energy_components("Sim Anneal", dr_state.step, state.energy, dr_state.temperature))
end

print_log_energy = Common.@callback 1000 function (state, dr_state)
    write(log_energy, print_energy_components("Sim Anneal", dr_state.step, state.energy, dr_state.temperature)*"\n")
    flush(log_energy)
end

print_struct_debug = Common.@callback 5000 function (state, dr_state)
    Print.as_pdb(debug_file, state, metadata)
    flush(debug_file)
end

# -------------------------------- Driver -------------------------------------
n_steps = 25_000
i_temp  = 1000.0

function adjust_temperature(step::Int64)
    return -(i_temp/n_steps) * step + i_temp
end

smpl_sa_driver = Drivers.MonteCarlo.DriverConfig(
    sampler = inner_cycle_sampler,
    evaluator = inner_cycle_evaluator,
    temperature = adjust_temperature,
    n_steps = n_steps,
    callbacks = [print_status_sa, print_log_energy, print_struct_debug])