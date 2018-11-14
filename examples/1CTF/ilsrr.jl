using ProtoSyn
using Printf

# ------------------
# SET CONFIGURATION:
# ------------------
const n_outer_steps            = 10     # Number of inner cycles to perform;
const n_inner_steps            = 2500   # Number of steps on each inner cycle;
const print_sts_every_ic       = 100    # Print system status every n steps;
const print_str_every_ic       = 0      # Print structure to files every n steps;
const min_step_size            = 0.01   # Minimum step_size allowed;
const updt_step_size_every     = 0      # Update step_size every n steps;

const dihedral_p_mut           = 0.025
dihedral_step_size             = π/2
const dihedral_p_mut_prt       = 0.25   # For perturbator
const dihedral_step_size_prt   = 0.01   # For perturbator

const crankshaft_p_mut         = 0.007
crankshaft_step_size           = π/2
const crankshaft_p_mut_prt     = 0.07   # For perturbator
const crankshaft_step_size_prt = 0.01   # For perturbator

const temperature              = 300.0  # Temperature for Metropolis acceptance in inner cycles (0.0 performs a Random Search where only better structures are accepted);
const acceptance_ratio         = 0.2    # Target acceptance ration of new structures. Step_size will be updated to try and achieve this value;
const ar_buffer_zone           = 0.01   # Buffer zone limit around the target acceptance_ratio. Step_size will not be updated if it's currently inside the buffer zone;
const hb_acceptance_criteria   = 1000.0 # Temperature for Metropolis acceptance in outer cycles. The bigger, the easier to accept worse configurations as new homebases;

const min_n_steps              = 25     # Number of steps for Steepest Descent minimization (every step of inner cycle);
const init_min_n_steps         = 10000  # Number of steps for Steepest Descent minimization (initial long minimization);
const print_sts_every_min      = 250    # Print status every n steps in Steepest Descent minimization;

const ss                       = "CEEEEEEECCCCHHHHHHHHHHHHCCCHHHHHHHHHCCCEEEEEEECHHHHHHHHHHHHHHCCEEEEC"
const input_pdb                = "data/1ctf_no_sc.pdb"
const input_contact_map        = "data/gremlin_contact_map.txt"
const input_amber_json         = "data/1ctf_amber_top.json"
const log_destination          = stdout
const xyz_destination          = open("out/trajectory.pdb", "w")
const best_destination         = open("out/best_trajectory.pdb", "w")
const out_destination          = open("out/data.out", "w")
const outer_best_destination   = open("out/outer_best_destination.pdb", "w")
const status_destination       = [log_destination, out_destination]

# ------------------
# DEFINE THE CALLBACKS:
# ------------------
include("ilsrr_callbacks.jl")


# ------------------
# INITIAL STEPS:
# ------------------
# Load state
state, metadata    = Common.load_from_pdb(input_pdb)

#Fix proline
Common.fix_proline!(state, metadata.dihedrals)

# Apply secondary structure
Common.apply_ss!(state, metadata, ss)

# Create blocks -> Will rotate side-chains, if present
nb_dihedrals       = filter(x -> x.residue.ss == Common.SS.COIL, metadata.dihedrals)
nb_phi_dihedrals   = filter(x -> x.dtype == Common.DIHEDRAL.phi, nb_dihedrals) # For crankshaft movements

# ------------------
# DEFINE THE DRIVERS:
# ------------------
include("ilsrr_drivers.jl")


# ------------------
# MAIN BODY:
# ------------------
# Define starting inner_best, outer_best and homebase *after minimization*
init_sd_driver = Drivers.SteepestDescent.SteepestDescentDriver(my_sd_evaluator!, n_steps = init_min_n_steps)
Drivers.SteepestDescent.run!(state, init_sd_driver, print_status_sd)
my_evaluator!(state, false)
inner_best = deepcopy(state)
outer_best = deepcopy(state)
homebase   = deepcopy(state)
Print.as_pdb(best_destination, state, metadata, title = "Inner best", step = 1)

# Run outer cycle
for outer_step in 1:n_outer_steps
    global state, inner_best, outer_best, homebase
    
    # Run inner cycle
    Print.status(@sprintf("\n(%5s) %12s \n%s\n", "ILSRR", @sprintf("Step: %4d", outer_step), "-"^110), log_destination)
    cb_status(0, state, mc_driver, 1.0, Dict("d" => 0, "c" => 0))
    Drivers.MonteCarlo.run!(state, mc_driver, print_status_mc, save_inner_best, adjust_step_size, print_structure, adjust_temperature)
    state = deepcopy(inner_best) # Output of the inner_cycle is always the inner_best
    
    # Save outer best
    if state.energy.eTotal < outer_best.energy.eTotal
        outer_best = deepcopy(state)
        Print.as_pdb(outer_best_destination, state, metadata, title = "Outer best", step = outer_step)
        flush(outer_best_destination)
        Print.status(@sprintf("(%5s) %12s | ⚡E: %10.3e | ✔️ Accepted\n", "ILSRR", "Inner Best", state.energy.eTotal), status_destination)
    else
        Print.status(@sprintf("(%5s) %12s | ⚡E: %10.3e | ✖️️ Not accepted\n", "ILSRR", "Inner Best", state.energy.eTotal), status_destination)
    end

    if (state.energy.eTotal < homebase.energy.eTotal) || (rand() < exp(-(state.energy.eTotal - homebase.energy.eTotal) / hb_acceptance_criteria))
        # Save new homebase (Pseudo Metropolis -> A worse configuration can still be a new homebase with a probability controled by `hb_acceptance_criteria`)
        homebase = deepcopy(state)
        Print.status("(ILSRR) New homebase defined ✔️\n", status_destination)
    else
        # Otherwise, recover the previous homebase for perturbation
        state = deepcopy(homebase)
    end

    # Perturb
    if outer_step != n_outer_steps
        mov = my_perturbator!(state)
        my_evaluator!(state, false)
        inner_best = deepcopy(state)
        s = @sprintf("(%5s) %12s | ⚡E: %10.3e | Dihedral: %3d | Crankshaft: %3d |\n%s\n", "ILSRR", "Perturbator", state.energy.eTotal, mov["d"], mov["c"], "-"^110)
        Print.status(s, status_destination)
        mc_driver.temperature = temperature
    end
end

close(xyz_destination)
close(best_destination)
close(outer_best_destination)