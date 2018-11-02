using ProtoSyn
using Printf

# ------------------
# SET CONFIGURATION:
# ------------------
const n_outer_steps            = 2
const n_inner_steps            = 20000
const print_every_inner_cycle  = 10

const dihedral_p_mut           = 0.02
dihedral_step_size             = π/2
const dihedral_p_mut_prt       = 1.0   # For perturbator
const dihedral_step_size_prt   = 0.01  # For perturbator

const crankshaft_p_mut         = 0.005
crankshaft_step_size           = 1.0
const crankshaft_p_mut_prt     = 1.0   # For perturbator
const crankshaft_step_size_prt = 0.01  # For perturbator
const crankshaft_min_n_steps   = 25    # For Steepest Descent mini-minimization

const temperature              = 400.15
const acceptance_ratio         = 0.2
const hb_acceptance_criteria   = 500.0 # The bigger, the easier to accept "worse" configurations as new homebases
const prt_min_n_steps          = 1000  # For Steepest Descent minimization after perturbator
const print_every_prt_minim    = 100   # For Steepest Descent minimization after perturbator

const ss                       = "CEEEEEEECCCCHHHHHHHHHHHHCCCHHHHHHHHHCCCEEEEEEECHHHHHHHHHHHHHHCCEEEEC"
const input_pdb                = "data/1ctf.pdb"
const input_contact_map        = "data/gremlin_contact_map.txt"
const input_amber_json         = "data/1ctf_amber_top.json"
const log_destination          = stdout
const xyz_destination          = open("out/trajectory.pdb", "w")
const best_destination         = open("out/best_trajectory.pdb", "w")
const out_destination          = open("out/data.out", "w")


# ------------------
# INITIAL STEPS:
# ------------------
# Load state
state, metadata    = Common.load_from_pdb(input_pdb)

#Fix proline
Common.fix_proline!(state, metadata.dihedrals)

# Apply secondary structure
Common.apply_ss!(state, metadata, ss)

# Create blocks -> Will rotate side-chains
nb_dihedrals       = filter(x -> x.residue.ss == Common.SS.COIL, metadata.dihedrals)
nb_phi_dihedrals   = filter(x -> x.dtype == Common.DIHEDRAL.phi, nb_dihedrals) # For crankshaft movements

# Define mutators
dihedral_mutator   = Mutators.Dihedral.DihedralMutator(nb_dihedrals, () -> (rand() * 2 - 1 * dihedral_mutator.step_size), dihedral_p_mut, dihedral_step_size)
crankshaft_mutator = Mutators.Crankshaft.CrankshaftMutator(nb_phi_dihedrals, () -> (rand() * 2 - 1 * crankshaft_mutator.step_size), crankshaft_p_mut, crankshaft_step_size)

# Define the Monte Carlo and Steepest Descent evaluator
topology           = Forcefield.Amber.load_from_json(input_amber_json)
contact_pairs      = Forcefield.Other.load_contact_maps_from_file(input_contact_map, metadata)
k = 100
function my_evaluator!(st::Common.State, do_forces::Bool)
    energy  = Forcefield.Amber.evaluate!(topology, st, cut_off=1.2, do_forces=do_forces)
    # println("AMBER: $energy")
    energy += k * Forcefield.Other.calc_eContact!(st, contact_pairs)
    # println("CONTACTS: $(k * st.energy.comp["eContact"])")
    # println("TOTAL: $energy")
    st.energy.eTotal = energy
    return energy
end

# Define the Monte Carlo sampler
sd_driver = Drivers.SteepestDescent.SteepestDescentDriver(my_evaluator!, n_steps = crankshaft_min_n_steps)
function my_sampler!(st::Common.State)
    dm = Mutators.Dihedral.run!(st, dihedral_mutator)
    cm = Mutators.Crankshaft.run!(st, crankshaft_mutator)
    # cm = 0

    # Crankshaft movements require a minimization step
    # if cm > 0 
    #     Drivers.SteepestDescent.run!(st, sd_driver)
    # end

    return Dict("Dihedral" => dm, "Crankshaft" => cm)
end

# Define the callbacks
# 1. Print status
dr_type = Union{Drivers.MonteCarlo.MonteCarloDriver, Drivers.SteepestDescent.SteepestDescentDriver}
function cb_status(step::Int64, st::Common.State, dr::Drivers.MonteCarlo.MonteCarloDriver, args...)
    write(log_destination, @sprintf "(%5s) %12d | ⚡E: %10.3e | Dihedral: %2d ▶️ %8.2e | Crankshaft: %2d ▶️ %8.2e | AR: %.4f | T: %4.2f\n" "MC" step st.energy.eTotal args[2]["Dihedral"] dihedral_mutator.step_size args[2]["Crankshaft"] crankshaft_mutator.step_size args[1] dr.temperature)
    write(out_destination, @sprintf "(%5s) %12d | ⚡E: %10.3e | Dihedral: %2d ▶️ %8.2e | Crankshaft: %2d ▶️ %8.2e | AR: %.4f | T: %4.2f\n" "MC" step st.energy.eTotal args[2]["Dihedral"] dihedral_mutator.step_size args[2]["Crankshaft"] crankshaft_mutator.step_size args[1] dr.temperature)
    flush(out_destination)
end
print_status_mc = Common.CallbackObject(print_every_inner_cycle, cb_status) 
print_status_sd = @Common.callback print_every_prt_minim function cb_status(step::Int64, st::Common.State, dr::Drivers.SteepestDescent.SteepestDescentDriver, args...)
    write(log_destination, @sprintf "(%5s) %12d | ⚡E: %10.3e | Max Force: %10.3e | Gamma: %10.3e\n" "SD" step st.energy.eTotal args[1] args[2])
end

# 2. Check if the structure is better than the current inner_best and save it
save_inner_best = @Common.callback 1 function cb_save(step::Int64, st::Common.State, dr::dr_type, args...)
    global inner_best
    if st.energy.eTotal < inner_best.energy.eTotal
        inner_best = deepcopy(st)
    end
end

# 3. Adjust the step_size so that, on average, the acceptance_ration is as defined initially
history = 0
adjust_step_size = @Common.callback 1 function cb_adjust_step_size(step::Int64, st::Common.State, dr::Drivers.MonteCarlo.MonteCarloDriver, ac, args...)
    global history
    d = ac > acceptance_ratio ? 1.0005 : 0.9995
    dihedral_mutator.step_size = max(1e-5, min(dihedral_mutator.step_size * d, π))
    crankshaft_mutator.step_size = max(1e-5, min(crankshaft_mutator.step_size * d, π))
    if dihedral_mutator.step_size == 1e-5
        history += 1
    else
        history = 0
    end
    if history > 10
        exit(1)
    end
end

# 4. Print current structure to a PDB file
print_structure = @Common.callback 100 function cb_print(step::Int64, st::Common.State, dr::dr_type, args...)
    Print.as_pdb(xyz_destination, st, metadata, step = step)
end

# Define the Monte Carlo Driver
mc_driver = Drivers.MonteCarlo.MonteCarloDriver(my_sampler!, my_evaluator!, temperature=temperature, n_steps=n_inner_steps)

# Define starting inner_best, outer_best and homebase
my_evaluator!(state, false)
inner_best = deepcopy(state)
outer_best = deepcopy(state)
homebase   = deepcopy(state)

# Define the perturbator
dihedral_perturbator = Mutators.Dihedral.DihedralMutator(nb_dihedrals, () -> (rand() * 2 - 1 * dihedral_mutator.step_size), dihedral_p_mut_prt, dihedral_step_size_prt)
crankshaft_perturbator = Mutators.Crankshaft.CrankshaftMutator(nb_phi_dihedrals, () -> (rand() * 2 - 1 * crankshaft_mutator.step_size), crankshaft_p_mut_prt, crankshaft_step_size_prt)
sd_prt_driver = Drivers.SteepestDescent.SteepestDescentDriver(my_evaluator!, n_steps = prt_min_n_steps)


# ------------------
# MAIN BODY:
# ------------------
for outer_step in 1:n_outer_steps
    global state, inner_best, outer_best, homebase

    #Inner cycle
    write(log_destination, @sprintf("\n(%5s) %12s \n%s\n", "ILSRR", @sprintf("Step: %4d", outer_step), "-"^110))
    Drivers.MonteCarlo.run!(state, mc_driver, print_status_mc, save_inner_best, adjust_step_size, print_structure)
    state = deepcopy(inner_best) # Output of the inner_cycle is always the inner_best
    println("INNER BEST: $(inner_best.energy.eTotal)")
    exit(1)
    
    # Save outer best
    if state.energy.eTotal < outer_best.energy.eTotal
        outer_best = deepcopy(state)
        Print.as_pdb(best_destination, state, metadata)
        write(log_destination, @sprintf "(%5s) %12s | ⚡E: %10.3e | ✔️ Accepted\n" "ILSRR" "Inner Best" state.energy.eTotal)
    else
        write(log_destination, @sprintf "(%5s) %12s | ⚡E: %10.3e | ✖️ Not accepted\n" "ILSRR" "Inner Best" state.energy.eTotal)
    end

    # Save new homebase (Pseudo Metropolis -> A worse configuration can still be a new homebase with a probability controled by `hb_acceptance_criteria`)
    if (state.energy.eTotal < homebase.energy.eTotal) || (rand() < exp(-(state.energy.eTotal - homebase.energy.eTotal) / hb_acceptance_criteria))
        homebase = deepcopy(state)
    end

    # Perturb
    if outer_step != n_outer_steps
        dm = Mutators.Dihedral.run!(state, dihedral_perturbator)
        cm = Mutators.Crankshaft.run!(state, crankshaft_perturbator)
        Drivers.SteepestDescent.run!(state, sd_prt_driver, print_status_sd)
        inner_best = deepcopy(state)
        write(log_destination, @sprintf "(%5s) %12s | ⚡E: %10.3e | Dihedral: %3d | Crankshaft: %3d |\n%s\n" "ILSRR" "Perturbator" state.energy.eTotal dm cm "-"^96)
    end
end

close(xyz_destination)
close(best_destination)