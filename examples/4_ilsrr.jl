using ProtoSyn
using Printf

#= -------------------------------------------------
4. ILSRR Example
(Iterated Local Search with Random Restarts)

# Algorithm explanation:
The ILSRR algorithm searchs the local conformational neighbourhood using
an "inner driver", such as Monte Carlo Simulated Annealing or even
Steepest Descent, and then "perturbates" the system in an "outer cycle",
in order to "kick" a solution out of the local optimum, essentially
overcomming energy barriers.

# New ProtoSyn concepts:
- Including reused code
- ILSRR & Compound drivers
-------------------------------------------------=#

# Include
#=
* Note: Whole files of code (containing configuration parameters,
callbacks, drivers, etc) can be included in another script, essentially
copying all content as if it was literally written here. Altough useful
in reusing code,This means all variable names are mainteined and caution
should be help not to overwrite or redefine essential pieces of code.
=#
include("3_main_sim_anneal.jl")

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
    @printf("%15s: %5d: %10.3f\n", "ILSRR", dr_state.step, dr_state.best_state.energy.total)
end


# ILSRR & Coumpound Drivers
#=
* Note: In this example, ILSRR Driver is a compound driver:
This means that one or more of its parameters are whole drivers. In this
example, the ILSRR Driver accepts a Monte Carlo Simulated Annealing driver
as its `inner_driver_config`. Both drivers can have diferent parameters,
such as the temperature or number of steps, aswell as different callbacks.
=#
ilsrr_driver = Drivers.ILSRR.DriverConfig(
    inner_driver_config = sa_driver,
    perturbator         = custom_perturbator,
    temperature         = 8.0,
    n_steps             = 1_000,
    stall_limit         = 20,
    callbacks           = [print_status_ilsrr]
)

if ""!=PROGRAM_FILE && realpath(@__FILE__) == realpath(PROGRAM_FILE)
    println("ILSRR Driver Example:")
    @time Drivers.run!(state, ilsrr_driver)
end