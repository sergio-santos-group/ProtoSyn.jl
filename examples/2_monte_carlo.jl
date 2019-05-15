using ProtoSyn
using Printf

#= -------------------------------------------------
2. Monte Carlo Example

# Algorithm explanation:
The Monte Carlo algorithm samples the conformational space, mutating
the structure or nature of the molecules in the simulation. Configurations
are visited more or less often based on the total energy of the system.
This energy can be a sum of different energy components.

# New ProtoSyn concepts:
- Energy components & Evaluators
- Samplers
- Callbacks
- Drivers
- Monte Carlo driver
-------------------------------------------------=#

# Configuration
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


#Energy Components & Evaluators
#=
* Note: In this example, the energy will be dictated by several components:
Amber forcefield, coarse grain approximations of solvation energy and
hydrogen bonding network, and distance restrains defined by long-range 
contact predictions. These components will be assembled in an Evaluator. By
default, this aggregation is cumulative: all components will be added to define
the total energy of the system. Custom aggregation functions may be defined.

An energy component is a static topology that defines for each atom how it
should behave and interact with its neighbours (as an example, the Amber topology
can dictate all ditances and force constants for bonds between consecutive atoms,
among others).
=#

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
    λ = 100.0,
    threshold = 0.0,
    min_distance = 0.8)

# Evaluators: An evaluator aggregates all defined components
custom_evaluator = Forcefield.Evaluator(
    components = [amber_top, solv_pairs, hb_network, contact_restraints])

# Sampler
#=
* Note: A Sampler is an aggregator of Mutators. By default, all mutators
defined are applied, one by one, but custom `apply!` function may be defined.
Additionally, a Sampler may receive an optional `tune!` function, which is
called in a specific point in the Driver algorithm in order to adjust the
mutators parameters (such as the step_size, as in this example).
=#

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
        mutator.step_size *= δ
    end
end
    
custom_sampler = Mutators.Sampler(
    mutators = [dihedral_mutator, crankshaft_mutator],
    tune! = adjust_step_size
)


#Callbacks
#=
* Note: Callbacks are functions with a peridiocity that are called inside
Drivers in order to comunicate with the user. These functions are called
at the end of each internal cycle of the Driver. In general, callbacks should
be unidirection, that is, they should comunicate the internal workings of
the algorithms outside, not changing any of the inner workings of the driver. 

Callbacks, by default, receive 2 arguments:
1. The current state (including coordinates, forces, energy, etc)
2. The driver state (including list of mutators, evaluator, step sizes, etc)

Multiple callbacks can be passed to Drivers. In this example they will comunicate
the current state of the system to the user every 100 steps while also printing
the current structure to a file, in a PDB format.
=#
print_status = Common.@callback 100 function (state, dr_state)
    @printf("%15s: %5d: %10.3f\n", "Monte Carlo", dr_state.step, state.energy.total)
end

print_structure = Common.@callback 100 function (state, dr_state)
    Print.as_pdb(output_file, state, metadata)
end

# Main Drivers
#=
* Note: Altough users are free to define their own custom algorithms,
ProtoSyn has several pre-defined algorithms. In this example, the Monte Carlo
will be used to sample the structure of a protein. Different drivers require
different input parameters, but overall a Driver usually has 4 necessary components:

1. A sampler: a motor to drive the state from one configuration to another. Is usually
comprised of a mix of mutators.
2. An evaluator: evaluates the energy of the system in the current configuration,
dictating the direction the Driver takes. Is usually a combination of different
energy components.
3. Configuration variables: configure various parameters of the driver itself, such as
the temperature of the system or number os steps, for example. Altough optional,
the default values for these parameters make it so an "empty" driver will usually not
produce any work, per se. 
4. Callbacks: comunicate the internal working of the driver to the outside or adjust
parameters during runtime. These are usually optional.

=#
mc_driver = Drivers.MonteCarlo.DriverConfig(
    sampler = custom_sampler,
    evaluator = custom_evaluator,
    temperature = 7.0,
    n_steps = 10_000,
    callbacks = [print_status, print_structure])

if ""!=PROGRAM_FILE && realpath(@__FILE__) == realpath(PROGRAM_FILE)
    const output_file = open("2_monte_carlo.pdb", "w")

    println("Monte Carlo example:")
    @time Drivers.run!(state, mc_driver)
    
    close(output_file)
end