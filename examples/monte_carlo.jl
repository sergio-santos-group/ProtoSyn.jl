# In this example, we will explore how to build a Monte Carlo simulation from
# scratch, using the Mutators and Drivers made available by ProtoSyn. As an
# example, we will try to predict the folded state of a sequence of aminoacids
# belonging to the 2A3D peptide from PDB.

using ProtoSyn
using ProtoSyn.Builder
using ProtoSyn.Peptides
using Printf

# 1) Load the residue library and build the peptide in a linear conformation
T        = Float64
res_lib  = Peptides.grammar(T)
sequence = seq"MGSWAEFKQRLAAIKTRLQALGGSEAELAAFEKEIAAFESELQAYKGKGNPEVEALRKEAAAIRDELQAYRHN"
begin
    pose     = Peptides.build(res_lib, sequence);

    # Secondary structure prediction increases the degrees of freedom or a
    # problem, by having more dihedral angles to predict. For this example, we will
    # use pre-calculated predictions for the secondary structure (from servers such
    # as RaptorX, among others), and apply the secondary structure to the respective
    # selections of residues.
    Peptides.setss!(pose, SecondaryStructure[:helix], rid"1:20");
    Peptides.setss!(pose, SecondaryStructure[:helix], rid"27:44");
    Peptides.setss!(pose, SecondaryStructure[:helix], rid"50:end");

    # For this stage, we won't be optimizing sidechain packaging. In order to
    # decrease the likelihood of steric clashes, we can't remove the sidechains
    # Peptides.remove_sidechains!(pose)
end

# 2) Load default energy function
energy_function = ProtoSyn.Common.get_default_energy_function()
energy_function.components[Peptides.Calculators.Caterpillar.solvation_energy] = 0.02
# energy_function = ProtoSyn.Calculators.EnergyFunction(Dict(
#     ProtoSyn.Calculators.TorchANI.torchani_model => T(1.0)), Int16(2))

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
p_mut              = 1/(length(selection(pose, gather = true)))
dihedral_mutator   = ProtoSyn.Mutators.DihedralMutator(
    randn, p_mut, 1.0, selection)

# For the CrankshaftMutator, only the CA atoms should be considered as rotation
# points. Furthermore, Proline residues should not be considered, given the
# ring-like backbone.
selection          = (rid"21:26" | rid"45:49") & an"CA" & !rn"PRO"
n                  = count(selection(pose))
p_mut              = 2/(n*(n-1))
crankshaft_mutator = ProtoSyn.Mutators.CrankshaftMutator(
    randn, p_mut, 2.0, selection, an"^C$|^O$"r)

# The two defined Mutators can now be combined in a CompoundDriver. This is
# simply an auxiliary object that iterates and calls any Function, Mutator or
# Driver in its body.
compound_driver = ProtoSyn.Drivers.CompoundDriver(
    [dihedral_mutator, crankshaft_mutator])

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
# quenching thermostat, starting at T = 0.1. Furthermore, optionally, a driver
# can receive a Callback function. This is called every step of the simulation
# and returns information to the user, such as the status of the pose or of the
# driver itself. For this example, we will use that information to save a
# trajectory of the simulation and print the current status to the user every 25
# steps. The simulation should run for 1000 steps.

function callback_function(pose::Pose, driver_state::ProtoSyn.Drivers.DriverState)
    @printf("STEP %-5d | E= %-10.4f | AR= %-5.1f%%  | T= %-7.4f\n", driver_state.step, pose.state.e[:Total], (driver_state.acceptance_count/driver_state.step)*100, driver_state.temperature)
    driver_state.step == 0 && return
    ProtoSyn.append(pose, "monte_carlo.pdb", model = driver_state.step)
end

callback    = ProtoSyn.Drivers.Callback(callback_function, 1000)
n_steps     = 100_000
monte_carlo = ProtoSyn.Drivers.MonteCarlo(
    energy_function,
    dihedral_mutator,
    callback,
    n_steps,
    ProtoSyn.Drivers.get_linear_quench(0.075, n_steps))
    # ProtoSyn.Drivers.get_constant_temperature(0.008))

# 5) Launch a simulation replica
begin
    ProtoSyn.write(pose, "monte_carlo.pdb")
    monte_carlo(pose)
end