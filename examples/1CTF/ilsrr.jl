using ProtoSyn
using Printf

# ------------------
# SET CONFIGURATION:
# ------------------
n_outer_steps            = 2
n_inner_steps            = 10000
print_every_inner_cycle  = 100

dihedral_p_mut           = 0.04
dihedral_step_size       = 0.001
dihedral_p_mut_prt       = 1.0   # (For perturbator
dihedral_step_size_prt   = 0.01  # For perturbator

crankshaft_p_mut         = 0.0005
crankshaft_step_size     = 0.01
crankshaft_p_mut_prt     = 1.0   # For perturbator
crankshaft_step_size_prt = 0.01  # For perturbator

temperature              = 0.2
acceptance_ratio         = 0.2
hb_acceptance_criteria   = 500.0 # The bigger, the easier to accept "worse" configurations as new homebases

ss                       = "CEEEEEEECCCCHHHHHHHHHHHHCCCHHHHHHHHHCCCEEEEEEECHHHHHHHHHHHHHHCCEEEEC"
input_pdb                = "data/1ctf.pdb"
input_mc_json            = "data/1ctf_mc_top.json"
input_amber_json         = "data/1ctf_amber_top.json"
log_destination          = stdout
xyz_destination          = open("out/trajectory.pdb", "w")


# ------------------
# INITIAL STEPS:
# ------------------
# Load state
state = Common.load_from_pdb(input_pdb)
state.energy = Forcefield.Amber.Energy()

#Fix proline
mc_topology = Aux.read_JSON(input_mc_json)
dihedrals, residues = Common.load_topology(mc_topology)
Common.fix_proline!(state, dihedrals)

# Apply secondary structure
Common.apply_ss!(state, dihedrals, ss)

# Create blocks -> Will rotate side-chains
nb_dihedrals = filter(x -> Int(x.residue.ss) < 1, dihedrals)

# Define mutators
dihedral_mutator = Mutators.Dihedral.DihedralMutator(nb_dihedrals, () -> (rand() * 2 - 1 * dihedral_mutator.step_size), dihedral_p_mut, dihedral_step_size)
crankshaft_mutator = Mutators.Crankshaft.CrankshaftMutator(nb_dihedrals, () -> (rand() * 2 - 1 * crankshaft_mutator.step_size), crankshaft_p_mut, crankshaft_step_size)

# Define the Monte Carlo sampler
function my_sampler!(st::Common.State)
    Mutators.Dihedral.run!(st, dihedral_mutator)
    # Mutators.Crankshaft.run!(st, crankshaft_mutator)
end

# Define the Monte Carlo evaluator
topology = Forcefield.Amber.load_from_json(input_amber_json)
function my_evaluator!(st::Common.State, do_forces::Bool)
    energy = Forcefield.Amber.evaluate!(topology, st, cut_off=1.2, do_forces=do_forces)
    return energy
end

# Define the callbacks
# 1. Print status
# 2. Check if the structure is better than the current inner_best and save it
# 3. Adjust the step_size so that, on average, the acceptance_ration is as defined initially
print_status = @Common.callback print_every_inner_cycle function cb_status(step::Int64, st::Common.State, dr::Drivers.MonteCarlo.MonteCarloDriver, args...)
    write(log_destination, @sprintf "(%5s) Step: %4d | Energy: %.4ef | D step: %.3e | CS step: %.3e | AR: %.4f\n" "MC" step state.energy.eTotal dihedral_mutator.step_size crankshaft_mutator.step_size args[1])
end
save_inner_best = @Common.callback 1 function cb_save(step::Int64, st::Common.State, dr::Drivers.MonteCarlo.MonteCarloDriver, args...)
    global inner_best
    if st.energy.eTotal < inner_best.energy.eTotal
        inner_best = deepcopy(st)
    end
end
adjust_step_size = @Common.callback 1 function cb_adjust_step_size(step::Int64, st::Common.State, dr::Drivers.MonteCarlo.MonteCarloDriver, args...)
    dcf = 0.1 # Dihedral change factor
    ccf = 0.1 # Crankshaft change factor
    if args[1] > acceptance_ratio
        dihedral_mutator.step_size * (1.0 + dcf) < π ? dihedral_mutator.step_size *= (1.0 + dcf) : dihedral_mutator.step_size = π
        crankshaft_mutator.step_size * (1.0 + ccf) < π ? crankshaft_mutator.step_size *= (1.0 + ccf) : crankshaft_mutator.step_size = π
    else
        dihedral_mutator.step_size * (1.0 - dcf) > 1e-5 ? dihedral_mutator.step_size *= (1.0 - dcf) : dihedral_mutator.step_size = 1e-5
        crankshaft_mutator.step_size * (1.0 - ccf) > 1e-5 ? crankshaft_mutator.step_size *= (1.0 - ccf) : crankshaft_mutator.step_size = 1e-5
    end
end
print_structure = @Common.callback 100 function cb_print(step::Int64, st::Common.State, dr::Drivers.MonteCarlo.MonteCarloDriver, args...)
    Print.as_pdb(xyz_destination, st, step = step)
end

# Define the Monte Carlo Driver
mc_driver = Drivers.MonteCarlo.MonteCarloDriver(my_sampler!, my_evaluator!, temperature=temperature, n_steps=n_inner_steps)

# Define starting inner_best, outer_best and homebase
inner_best = deepcopy(state)
outer_best = deepcopy(state)
homebase   = deepcopy(state)

# Define the perturbator
dihedral_perturbator = Mutators.Dihedral.DihedralMutator(nb_dihedrals, () -> (rand() * 2 - 1 * dihedral_mutator.step_size), dihedral_p_mut_prt, dihedral_step_size_prt)
crankshaft_perturbator = Mutators.Crankshaft.CrankshaftMutator(nb_dihedrals, () -> (rand() * 2 - 1 * crankshaft_mutator.step_size), crankshaft_p_mut_prt, crankshaft_step_size_prt)


# ------------------
# MAIN BODY:
# ------------------
Print.as_pdb(xyz_destination, state)
Drivers.MonteCarlo.run!(state, mc_driver, print_status, save_inner_best, adjust_step_size, print_structure)
Print.as_pdb(xyz_destination, state)
exit(1)
for outer_step in 1:n_outer_steps
    global state, inner_best, outer_best, homebase

    #Inner cycle
    write(log_destination, @sprintf "\n(%5s) Step: %4d \n%s\n" "ILSRR" outer_step "-"^94)
    Drivers.MonteCarlo.run!(state, mc_driver, print_status, save_inner_best, adjust_step_size, print_structure)
    state = deepcopy(inner_best) # Output of the inner_cycle is always the inner_best

    # Save outer best
    println(state.xyz)
    println("$(state.energy.eTotal) < $(outer_best.energy.eTotal) ? $(state.energy.eTotal < outer_best.energy.eTotal)")
    if state.energy.eTotal < outer_best.energy.eTotal
        outer_best = deepcopy(state)
        Print.as_pdb(xyz_destination, state)
    end

    # Save new homebase (Pseudo Metropolis -> A worse configuration can still be a new homebase with a probability controled by `hb_acceptance_criteria`)
    if (state.energy.eTotal < homebase.energy.eTotal) || (rand() < exp(-(state.energy.eTotal - homebase.energy.eTotal) / hb_acceptance_criteria))
        homebase = deepcopy(state)
    end

    # Perturb
    Mutators.Dihedral.run!(state, dihedral_perturbator)
    Mutators.Crankshaft.run!(state, crankshaft_perturbator)
end