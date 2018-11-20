using ProtoSyn
using Printf

# SET CONFIGURATION:
# ------------------
include("ilsrr_config.jl")

# DEFINE THE CALLBACKS:
# ------------------
include("ilsrr_callbacks.jl")

# INITIAL STEPS:
# ------------------
# Load state and metadata
state, metadata    = Common.load_from_pdb(input_pdb)
#Fix proline
# Common.fix_proline!(state, metadata.dihedrals)
# Apply secondary structure
# Common.apply_ss!(state, metadata, ss)
metadata.ss = Common.infer_ss(metadata.dihedrals, ss)
# Create blocks -> Will rotate side-chains, if present
nb_dihedrals       = filter(x -> x.residue.ss == Common.SS.COIL, metadata.dihedrals)
nb_phi_dihedrals   = filter(x -> x.dtype == Common.DIHEDRAL.phi, nb_dihedrals) # For crankshaft movements

Print.as_pdb(xyz_destination, state, metadata, title = "Init", step = 1)

# DEFINE THE DRIVERS:
# ------------------
include("ilsrr_drivers.jl")
# MAIN BODY:
# ------------------
# Define starting inner_best, outer_best and homebase *after minimization*
init_sd_driver = Drivers.SteepestDescent.SteepestDescentDriver(my_sd_evaluator!, n_steps = init_min_n_steps)
Drivers.SteepestDescent.run!(state, init_sd_driver, print_status_sd, print_structure)
my_evaluator!(state, false)
inner_best = deepcopy(state)
outer_best = deepcopy(state)
homebase   = deepcopy(state)
Print.as_pdb(best_destination, state, metadata, title = "Inner best", step = 1)
exit(1)

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

    # Save / Recover homebase
    if (state.energy.eTotal < homebase.energy.eTotal) || (rand() < exp(-(state.energy.eTotal - homebase.energy.eTotal) / outer_temperature))
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
        mc_driver.temperature = inner_temperature
    end
end

close(xyz_destination)
close(best_destination)
close(outer_best_destination)