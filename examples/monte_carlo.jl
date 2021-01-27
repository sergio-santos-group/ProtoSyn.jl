# In this example, we will explore how to build a Monte Carlo simulation from
# scratch, using the Mutators and Drivers made available by ProtoSyn. As an
# example, we will try to loosely predict the folded state of a sequence of
# aminoacids belonging to the 2A3D peptide from PDB.

using ProtoSyn
using ProtoSyn.Builder
using ProtoSyn.Peptides
using Printf
using Random

# 1) Load the residue library and build the peptide in a linear conformation
T        = Float64
res_lib  = Peptides.grammar(T)
sequence = seq"MGSWAEFKQRLAAIKTRLQALGGSEAELAAFEKEIAAFESELQAYKGKGNPEVEALRKEAAAIRDELQAYRHN"

# We can assign a specific seed to replicate a previous run or perform a new
# random sampling simulation
begin
    seed = rand(1:1000000)
    # seed = 103148
    Random.seed!(seed)
    pose = Peptides.build(res_lib, sequence);

    # Secondary structure prediction increases the degrees of freedom of a problem,
    # by having more dihedral angles to predict. For this example, we will use
    # pre-calculated predictions for the secondary structure (from servers such as
    # RaptorX, among others), and apply the secondary structure
    # to the respective selections of residues.
    Peptides.setss!(pose, SecondaryStructure[:helix], rid"1:20");
    Peptides.setss!(pose, SecondaryStructure[:helix], rid"27:44");
    Peptides.setss!(pose, SecondaryStructure[:helix], rid"50:end");

    # For this stage, we won't be optimizing sidechain packaging. In order to 
    # decrease the likelihood of steric clashes, we can remove the sidechains
    Peptides.remove_sidechains!(pose)

    # We can also load the known biological structure to act as a comparison point
    # between the base truth and the simulated prediction
    base_truth = Peptides.load("../2a3d.pdb")
    Peptides.remove_sidechains!(base_truth)

    # 2) Load default energy function
    # After loading the default energy function we can change the relative weight of
    # each component. For this example, we will load an additional component: a
    # contact map. This map can be predicted using other software or online servers,
    # such as RaptorX, and gives the likelihood (in percentage) of two non
    # consecutive residues having a distance no greater than 8Å between them, based
    # on machine learning models trained on pre-existing curated databases of known
    # structures. This information is therefore more useful in 3D structure
    # prediction simulations (rather than design efforts) if homologue structures
    # exist and whose structure is known.
    energy_function = ProtoSyn.Common.default_energy_function()

    cm_restraints   = Peptides.Calculators.Restraints.load_contact_map("contact_map_example.txt")
    energy_function.components[ProtoSyn.Peptides.Calculators.Caterpillar.solvation_energy] = 0.001
    energy_function.components[cm_restraints] = 0.0001

    # 3) Define the Mutators
    # A Mutator is a piece of code responsible for performing a given conformational
    # change in the molecular structure. For this example, we will employ 2
    # different Mutators: a DihedralMutator and a CrankshaftMutator. In order to
    # minimize the degrees of freedom of the problem, and thus increase the speed of
    # the calculation, only a subset of residues will be the target for the
    # Mutators, more especifically, the loop areas (the pre-set secondary structure
    # section will be maintained). For the DihedralMutator, only the C and N atoms
    # of the backbone will be selected (these control the Phi and Psi angles).
    selection          = (rid"21:26" | rid"45:49") & an"C$|N$"r
    p_mut              = 1/(count(selection(pose)))
    dihedral_mutator   = ProtoSyn.Mutators.DihedralMutator(
        randn, p_mut, 0.5, selection)

    # For the CrankshaftMutator, only the CA atoms should be considered as rotation
    # points. Furthermore, Proline residues should not be considered, given the
    # ring-like backbone. Note: A crankshaft rotation can cause steric clashes,
    # especially on the last residue downstream of the selected region. In order to
    # minimize this, the step_size of this mutator should be relatively small.
    selection          = (rid"21:26" | rid"45:49") & an"CA" & !rn"PRO"
    n                  = count(selection(pose))
    p_mut              = 0.5/(n*(n-1))
    # p_mut = 0
    crankshaft_mutator = ProtoSyn.Mutators.CrankshaftMutator(
        randn, p_mut, 0.05, selection, !(an"^CA$|^N$|^C$|^H$|^O$"r))

    # The two defined Mutators can now be combined in a CompoundDriver. This is
    # simply an auxiliary object that iterates and calls any Function, Mutator or
    # Driver in its body.
    compound_driver = ProtoSyn.Drivers.CompoundDriver(
        [crankshaft_mutator, dihedral_mutator])

    # 4) Define the Monte Carlo Driver
    # A Driver is a more complex piece of code that is able to direct the flow of
    # the code, by making decisions based on the evaluation of a particular state of
    # the system. As an example, a Monte Carlo Driver performs a conformational
    # change, evaluates the resulting system and continues or returns to the
    # previous saved state based on the Metropolis Criterium (at a given
    # temperature). Besides the energy function and mutators, the Monte Carlo Driver
    # also needs a temperature function. These can be set by the user, but ProtoSyn
    # makes available a set of pre-defined "thermostats" that control the
    # temperature in some specific way. For this example, we will use a temperature
    # quenching thermostat, starting at T = 0.5. Furthermore, optionally, a driver
    # can receive a Callback function. This is called every step of the simulation
    # and returns information to the user, such as the status of the pose or of the
    # driver itself. For this example, we will use that information to save a
    # trajectory of the simulation and print the current status to the user every
    # 1000 steps, as well as to a .dat file. The simulation should run for 15 000
    # steps.

    function callback_function(pose::Pose, driver_state::ProtoSyn.Drivers.DriverState)
        s = @sprintf(" INNER STEP %-6d | E(total)= %-10.4f | AR= %-5.1f%%  | T= %-7.4f | E(ani)= %-10.4f | E(sol)= %-10.4f | E(cmp)= %-10.4f | E(cls)= %-10.4f\n", driver_state.step, pose.state.e[:Total], (driver_state.acceptance_count/driver_state.step)*100, driver_state.temperature, pose.state.e[:TorchANI_ML_Model], pose.state.e[:Caterpillar_Solvation], pose.state.e[:Contact_Map_Restraint], pose.state.e[:Clash_Restraint])
        print(s)
        open("monte_carlo.dat", "a") do io
            write(io, s)
        end
        driver_state.step == 0 && return
        ProtoSyn.append(pose, "monte_carlo.pdb")
    end

    callback    = ProtoSyn.Drivers.Callback(callback_function, 1000)
    n_steps     = 15_000
    monte_carlo = ProtoSyn.Drivers.MonteCarlo(
        energy_function,
        compound_driver,
        callback,
        n_steps,
        ProtoSyn.Drivers.get_linear_quench(0.01, n_steps))

    # 5) Define the steepest descent driver
    # In this example, we can relax the initial structure to conform to the
    # TorchANI model (all other components weights are set to zero or do not
    # contribute to forces calculation).
    sd_energy_function = copy(energy_function)
    sd_energy_function.components[cm_restraints] = 0.0

    # ProtoSyn makes available certain common default callbacks. In this case,
    # the "energy_step" callback simply prints the current step and total energy
    # of the system (every 500 steps, in this case).
    cb = ProtoSyn.Common.default_energy_step_callback(500)
    steepest_descent = ProtoSyn.Drivers.SteepestDescent(
        sd_energy_function, cb, 3000, 0.001, 0.1) # Note the energy function

    # 6) Launch a simulation replica

    # Print the base truth information
    println("Minimizing base truth:")
    steepest_descent(base_truth)
    base_truth_energy = energy_function(base_truth)
    open("monte_carlo.dat", "w") do io
        s = @sprintf("  BASE STEP %-6d | E(total)= %-10.4f | AR= %-5.1f%%  | T= %-7.4f | E(ani)= %-10.4f | E(sol)= %-10.4f | E(cmp)= %-10.4f | E(cls)= %-10.4f\n", 0, base_truth.state.e[:Total], 0.0, 0.0, base_truth.state.e[:TorchANI_ML_Model], base_truth.state.e[:Caterpillar_Solvation], base_truth.state.e[:Contact_Map_Restraint], base_truth.state.e[:Clash_Restraint])
        println(s)
        write(io, s)
    end

    # Start the simulation
    ProtoSyn.write(pose, "monte_carlo.pdb")
    steepest_descent(pose)
    monte_carlo(pose)

    println("Seed: $seed")
end