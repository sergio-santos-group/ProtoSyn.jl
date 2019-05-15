using ProtoSyn
using Printf

#= -------------------------------------------------
2. Simulated Annealing Example

# Algorithm explanation:
The Simulated Annealing algorithm searches the local conformational space,
mutating the structure or nature of the molecules in the simulation.
Configurations are accepted or not based on the Metropolis criterium.
However, differing from a simple Monte Carlo, the temperature of system
is gradually lowered until 0, so that, in theory, a single local optimum is
identified.
-------------------------------------------------=#

# Configuration
const input_pdb                = "data/1i2t_no_sc.pdb"
const input_amber_top          = "data/1i2t_amber_no_sc_top.json"
const hydrogen_bonding_library = "data/sc_hb_lib.json"
const contacts_topology        = "data/contact_map_raptorx_1i2t.txt"
const ss                       = "CHHHHHHHHHHHHHHHCCCCHHHHHHHHHHCCHHHHHHHHHCHHHHHHHHHHHHHHHHHHC"
const output_file              = open("3_sim_anneal.pdb", "w")

# System Set-up & Loading
state, metadata = Common.load_from_pdb(input_pdb)
Common.apply_ss!(state, metadata, ss)
nb_dhs          = filter(x -> x.residue.ss == Common.SS.COIL, metadata.dihedrals)
nb_phi_dhs      = filter(x -> x.dtype == Common.DIHEDRAL.phi, nb_dhs)

# Component A: Amber forcefield
amber_top = Forcefield.Amber.load_from_json(input_amber_top)

# Component B: Coarse-grain solvation energy
solv_pairs = Forcefield.CoarseGrain.compile_solv_pairs(
    metadata.dihedrals,
    λ = 1.0
)

# Component C: Coarse-grain hydrogen bonding network
hb_network = Forcefield.CoarseGrain.compile_hb_network(
    metadata.atoms,
    lib = Aux.read_JSON(hydrogen_bonding_library),
    λ = 1.0)

# Component D: Contacts distance restraints
contact_restraints = Forcefield.Restraints.load_distance_restraints_from_file(
    contacts_topology,
    metadata,
    λ = 2e3,
    threshold = 0.0,
    min_distance = 0.8)

# Evaluators: An evaluator aggregates all defined components
custom_evaluator = Forcefield.Evaluator(
    components = [amber_top, solv_pairs, hb_network, contact_restraints])

# Mutators
dihedral_mutator   = Mutators.Dihedral.MutatorConfig(
    dihedrals = nb_dhs,
    angle_sampler = () -> (randn() * dihedral_mutator.step_size),
    p_mut = 0.025,
    step_size = π/8)

crankshaft_mutator = Mutators.Crankshaft.MutatorConfig(
    dihedrals = nb_phi_dhs,
    angle_sampler = () -> (randn() * crankshaft_mutator.step_size),
    p_mut = 0.007,
    step_size = π/8)
    
function adjust_step_size(mutators, dr_state)
    ac_ratio = dr_state.ac_count / dr_state.step
    δ = ac_ratio > 0.2 ? 1.05 : 0.75
    for mutator in mutators
        mutator.step_size = max(0.001, min(mutator.step_size * δ, π))
    end
end
    
custom_sampler = Mutators.Sampler(
    mutators = [dihedral_mutator, crankshaft_mutator],
    tune! = adjust_step_size
)

# Callbacks
print_status = Common.@callback 100 function (st, dr_state)
    @printf("%15s: %5d: %10.3f (T %5.3f)\n", "Sim Annealing", dr_state.step, st.energy.total, dr_state.temperature)
end

print_structure = Common.@callback 100 function (st, dr_state)
    Print.as_pdb(output_file, st, metadata)
end

function adjust_temperature(step::Int64)
    return -0.007 * step + 7.0
end

sa_driver = Drivers.MonteCarlo.DriverConfig(
    sampler = custom_sampler,
    evaluator = custom_evaluator,
    temperature = adjust_temperature,
    n_steps = 1_000,
    callbacks = [print_status])

if ""!=PROGRAM_FILE && realpath(@__FILE__) == realpath(PROGRAM_FILE)
    println("Simulated Annealing example:")
    @time Drivers.run!(state, sa_driver)
end