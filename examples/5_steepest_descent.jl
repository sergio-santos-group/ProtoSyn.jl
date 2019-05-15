using ProtoSyn
using Printf

#= -------------------------------------------------
2. Steepest Descent Example

# Algorithm explanation:
Also known as Gradient Descent, Steepest Descent is a first order iterative
optimization algorithm, using the negative of the derivate of the energy
function to calculate the forces being felt by the system and applying them
in controled step sizes, essentially moving the system to a local minimum.

# New ProtoSyn concepts:
- SteepestDescent Driver
-------------------------------------------------=#

# Configuration
input_pdb                = "data/1i2t_no_sc_twist.pdb"
input_amber_top          = "data/1i2t_amber_no_sc_top.json"
ss                       = "CHHHHHHHHHHHHHHHCCCCHHHHHHHHHHCCHHHHHHHHHCHHHHHHHHHHHHHHHHHHC"


# System Set-up & Loading
state, metadata = Common.load_from_pdb(input_pdb)
Common.apply_ss!(state, metadata, ss)
nb_dhs          = filter(x -> x.residue.ss == Common.SS.COIL, metadata.dihedrals)
nb_phi_dhs      = filter(x -> x.dtype == Common.DIHEDRAL.phi, nb_dhs)

# Component A: Amber forcefield
amber_top = Forcefield.Amber.load_from_json(input_amber_top)

# Evaluators: An evaluator aggregates all defined components
custom_evaluator = Forcefield.Evaluator(
    components = [amber_top])

# Callbacks
print_status = Common.@callback 100 function (st::Common.State, dr_state)
    @printf("%15s %5d: %10.3f %10.3e\n", "Minimizing...", dr_state.step, st.energy.total, dr_state.step_size)
end

print_structure = Common.@callback 100 function (st::Common.State, dr_state)
    Print.as_pdb(output_file, st, metadata)
end

# Main Drivers
sd_driver = Drivers.SteepestDescent.DriverConfig(
    evaluator   = custom_evaluator,
    n_steps     = 2_000,
    nblist_freq = 10,
    f_tol       = 0.2,
    max_step    = 0.01,
    callbacks   = [print_status])

if ""!=PROGRAM_FILE && realpath(@__FILE__) == realpath(PROGRAM_FILE)
    const output_file = open("5_steepest_descent.pdb", "w")

    println("Steepest Descent example:")
    @printf("%15s %5s: %10s %10s\n", "", "Step", "Energy", "Step size")

    @time f_state = Drivers.run!(state, sd_driver)
    println(f_state)

    close(output_file)
end