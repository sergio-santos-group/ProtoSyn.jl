using ProtoSyn
using Printf

# ------------------ Aux functions -------------------------

function print_energy_components(driver::String, step::Int64, energy::Common.Energy, temp::Float64)::String
    components = [:amber, :contacts, :sol, :hb]
    s = @sprintf("%10s %6d %11.4e", driver, step, energy.total)
    for component in components
        s = join([s, @sprintf("%11.4e", energy.components[component])], " ")
    end
    s = join([s, @sprintf("| %6.3f", temp)], " ")
    return s
end

function print_energy_keys(driver::String, step::Int64, energy::Common.Energy)::String
    components = [:amber, :contacts, :sol, :hb]
    s = @sprintf("%10s %6d %11s", driver, step, "TOTAL")
    for key in components
        s = join([s, @sprintf("%11s", uppercase(string(key)))], " ")
    end
    s = join([s, " | TEMPERATURE"])
    return s
end

# ---------------------------------------------------------------

# System configuration
input_pdb                = "data/1i2t_no_sc.pdb"
input_amber_top          = "data/1i2t_amber_no_sc_top.json"
hydrogen_bonding_library = "data/sc_hb_lib.json"
contacts_topology        = "data/contact_map_raptorx_1i2t.txt"
ss                       = "CHHHHHHHHHHHHHHHCCCCHHHHHHHHHHCCHHHHHHHHHCHHHHHHHHHHHHHHHHHHC"

# System Set-up & Loading
state, metadata = Common.load_from_pdb(input_pdb)
Common.apply_ss!(state, metadata, ss)
nb_dhs          = filter(x -> x.residue.ss == Common.SS.COIL, metadata.dihedrals)
nb_phi_dhs      = filter(x -> x.dtype == Common.DIHEDRAL.phi, nb_dhs)

amber_top = Forcefield.Amber.load_from_json(input_amber_top)

# Component B: Coarse-grain solvation energy
solv_pairs = Forcefield.CoarseGrain.compile_solv_pairs(
    metadata.dihedrals,
    λ = 0.01)

# Component C: Coarse-grain hydrogen bonding network
hb_network = Forcefield.CoarseGrain.compile_hb_network(
    metadata.atoms,
    lib = Aux.read_JSON(hydrogen_bonding_library),
    λ = 2e4)

# Component D: Contacts distance restraints
contact_restraints = Forcefield.Restraints.load_distance_restraints_from_file(
    contacts_topology,
    metadata,
    λ = 20.0,
    threshold = 0.0,
    min_distance = 0.8)

# Evaluators: An evaluator aggregates all defined components
custom_evaluator = Forcefield.Evaluator(
    components = [amber_top, solv_pairs, hb_network, contact_restraints])

# ---------------------------------------------------
#       I N N E R    C Y C L E: Monte Carlo
# ...................................................

# Mutators
dihedral_mutator   = Mutators.Dihedral.MutatorConfig(
    dihedrals = nb_dhs,
    angle_sampler = () -> (randn() * dihedral_mutator.step_size),
    p_mut = 0.05,
    step_size = π/8)

crankshaft_mutator = Mutators.Crankshaft.MutatorConfig(
    dihedrals = nb_phi_dhs,
    angle_sampler = () -> (randn() * crankshaft_mutator.step_size),
    p_mut = 0.027,
    step_size = π/8)

function adjust_step_size(mutators, dr_state)
    ac_ratio = dr_state.ac_count / dr_state.step
    δ = ac_ratio > 0.2 ? 1.05 : 0.75
    for mutator in mutators
        mutator.step_size = max(0.001, min(mutator.step_size * δ, π))
    end
end

# Sampler
inner_cycle_sampler = Mutators.Sampler(
    mutators = [dihedral_mutator, crankshaft_mutator],
    tune! = adjust_step_size
)

# Inner cycle callbacks
print_status_sa = Common.@callback 5000 function (state, dr_state)
    println(print_energy_components("MonteCarlo", dr_state.step, state.energy, dr_state.temperature))

end

print_struct_sa = Common.@callback 10000 function (state, dr_state)
    Print.as_pdb(output_file, state, metadata)
    flush(output_file)
end

#            ~> Inner cycle driver <~
# `´``´``´``´``´``´``´``´``´``´``´``´``´``´``´``´``´``
n_steps = 100_000
i_temp = 70.0
function adjust_temperature(step::Int64)
    return -(i_temp/n_steps) * step + i_temp
end
sa_driver = Drivers.MonteCarlo.DriverConfig(
    sampler = inner_cycle_sampler,
    evaluator = custom_evaluator,
    temperature = adjust_temperature,
    n_steps = n_steps,
    callbacks = [print_status_sa, print_struct_sa])


# ---------------------------------------------------
#           O U T E R    C Y C L E: ILSRR
# ...................................................

# Perturbators
dihedral_perturbator   = Mutators.Dihedral.MutatorConfig(
    dihedrals = nb_dhs,
    angle_sampler = () -> (randn() * dihedral_mutator.step_size),
    p_mut = 0.25,
    step_size = π/2)

crankshaft_perturbator = Mutators.Crankshaft.MutatorConfig(
    dihedrals = nb_phi_dhs,
    angle_sampler = () -> (randn() * crankshaft_mutator.step_size),
    p_mut = 0.07,
    step_size = π/2)

# Perturbator Sampler (Perturbator aggregator)
custom_perturbator = Mutators.Sampler(
    mutators = [dihedral_perturbator, crankshaft_perturbator])

print_status_ilsrr = Common.@callback 1 function (state, dr_state)
    printstyled(print_energy_keys("ILSRR", dr_state.step, dr_state.best_state.energy)*"\n", color = :green)
    printstyled(print_energy_components("ILSRR", dr_state.step, dr_state.best_state.energy, dr_state.temperature)*"\n", color = :green)
end

print_outer_best = Common.@callback 1 function (state, dr_state)
    Print.as_pdb(outer_best_file, dr_state.best_state, metadata)
    flush(outer_best_file)
end

#            ~> Outer cycle driver <~
# `´``´``´``´``´``´``´``´``´``´``´``´``´``´``´``´``´``
ilsrr_driver = Drivers.ILSRR.DriverConfig(
    inner_driver_config = sa_driver,
    perturbator         = custom_perturbator,
    temperature         = 70.0,
    n_steps             = 10_000,
    stall_limit         = 20,
    callbacks           = [print_status_ilsrr, print_outer_best])



# Main Run <~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if ""!=PROGRAM_FILE && realpath(@__FILE__) == realpath(PROGRAM_FILE)
    const output_file = open("inner_cycle.pdb", "w")
    const outer_best_file = open("outer_cycle.pdb", "w")

    println("SAMPLER:")
    @time Drivers.run!(state, ilsrr_driver)
end