using ProtoSyn
using ProtoSyn
using ProtoSyn.Peptides
using Printf
using Random

T        = Float64
res_lib  = Peptides.grammar(T)
sequence = seq"MGSWAEFKQRLAAIKTRLQALGGSEAELAAFEKEIAAFESELQAYKGKGNPEVEALRKEAAAIRDELQAYRHN"

base_truth = Peptides.load("../2a3d.pdb")
Peptides.remove_sidechains!(base_truth)

begin
    seed = rand(1:1000000)
    # seed = 150223 (5.27Å if α(cm_restraints) = 0.01)
    Random.seed!(seed)
    pose = Peptides.build(res_lib, sequence);

    Peptides.setss!(pose, SecondaryStructure[:helix], rid"1:20");
    Peptides.setss!(pose, SecondaryStructure[:helix], rid"27:44");
    Peptides.setss!(pose, SecondaryStructure[:helix], rid"50:end");

    Peptides.remove_sidechains!(pose)
end

begin
    energy_function = ProtoSyn.Common.default_energy_function()
    energy_function.components[Peptides.Calculators.Caterpillar.solvation_energy] = 0.01

    cm_restraints   = Peptides.Calculators.Restraints.load_contact_map("contact_map_example.txt")
    energy_function.components[cm_restraints] = 0.001

    sd_energy_function = copy(energy_function)
    sd_energy_function.components[cm_restraints] = 0.0
end


# Define the inner monte carlo driver
begin
    selection              = (rid"21:26" | rid"45:51") & an"C$|N$"r
    p_mut                  = 1/(count(selection(pose)))
    inner_dihedral_mutator = ProtoSyn.Mutators.DihedralMutator(
        randn, p_mut, 1.0, selection)

    selection                = (rid"21:26" | rid"45:51") & an"CA" & !rn"PRO"
    n                        = count(selection(pose))
    p_mut                    = 0.5/(n*(n-1))
    inner_crankshaft_mutator = ProtoSyn.Mutators.CrankshaftMutator(
        randn, p_mut, 0.1, selection, !(an"^CA$|^N$|^C$|^H$|^O$"r))

    inner_compound_driver    = ProtoSyn.Drivers.CompoundDriver(
        [inner_crankshaft_mutator, inner_dihedral_mutator])

    function inner_callback_function(pose::Pose, driver_state::ProtoSyn.Drivers.DriverState)
        s = @sprintf(" INNER STEP %-6d | E(total)= %-10.4f | AR= %-5.1f%%  | T= %-7.4f | E(ani)= %-10.4f | E(sol)= %-10.4f | E(cmp)= %-10.4f | E(cls)= %-10.4f\n", driver_state.step, pose.state.e[:Total], (driver_state.acceptance_count/driver_state.step)*100, driver_state.temperature, pose.state.e[:TorchANI_ML_Model], pose.state.e[:Caterpillar_Solvation], pose.state.e[:Contact_Map_Restraint], pose.state.e[:Clash_Restraint])
        print(s)
        open("ilsrr.dat", "a") do io
            write(io, s)
        end
        driver_state.step == 0 && return
        ProtoSyn.append(pose, "ilsrr.pdb")
    end

    inner_callback    = ProtoSyn.Drivers.Callback(inner_callback_function, 1000)
    inner_n_steps     = 10_000
    inner_monte_carlo = ProtoSyn.Drivers.MonteCarlo(
        energy_function,
        inner_compound_driver,
        inner_callback,
        inner_n_steps,
        ProtoSyn.Drivers.get_linear_quench(0.3, inner_n_steps))
end


# Define the outer ILSRR driver
begin
    selection              = (rid"21:26" | rid"45:51") & an"C$|N$"r
    outer_dihedral_mutator = ProtoSyn.Mutators.DihedralMutator(
        randn, 0.8, 0.2, selection)

    selection                = (rid"21:26" | rid"45:51") & an"CA" & !rn"PRO"
    n                        = count(selection(pose))
    outer_crankshaft_mutator = ProtoSyn.Mutators.CrankshaftMutator(
        randn, 0.15, 0.05, selection, !(an"^CA$|^N$|^C$|^H$|^O$"r))

    cb = ProtoSyn.Common.default_energy_step_callback(500)
    steepest_descent         = ProtoSyn.Drivers.SteepestDescent(
        sd_energy_function, cb, 3000, 0.001, 0.1) # Note the energy function

    outer_compound_driver    = ProtoSyn.Drivers.CompoundDriver(
        [outer_crankshaft_mutator, outer_dihedral_mutator, steepest_descent])

    function outer_callback_function(pose::Pose, driver_state::ProtoSyn.Drivers.DriverState)
        s = @sprintf("OUTER STEP %-6d | E= %-10.4f | AR= %-5.1f%%  | T= %-7.4f\n", driver_state.step, pose.state.e[:Total], (driver_state.acceptance_count/driver_state.step)*100, driver_state.temperature)
        print(s)
        open("ilsrr.dat", "a") do io
            write(io, s)
        end
        driver_state.step == 0 && return
        ProtoSyn.append(pose, "ilsrr_best.pdb")
    end

    outer_callback    = ProtoSyn.Drivers.Callback(outer_callback_function, 1)
    outer_n_steps     = 100
    outer_ILSRR = ProtoSyn.Drivers.ILSRR(
        energy_function,
        outer_compound_driver,
        inner_monte_carlo,
        outer_callback,
        outer_n_steps,
        ProtoSyn.Drivers.get_linear_quench(0.1, outer_n_steps))
end

begin
    base_truth_energy = energy_function(base_truth)
    open("ilsrr.dat", "w") do io
        s = @sprintf(" BASE STEP %-6d | E(total)= %-10.4f | AR= %-5.1f%%  | T= %-7.4f | E(ani)= %-10.4f | E(sol)= %-10.4f | E(cmp)= %-10.4f\n", 0, base_truth.state.e[:Total], 0.0, 0.0, base_truth.state.e[:TorchANI_ML_Model], base_truth.state.e[:Caterpillar_Solvation], base_truth.state.e[:Contact_Map_Restraint])
        println(s)
        write(io, s)
    end
    sync!(pose)
    ProtoSyn.write(pose, "ilsrr.pdb")
    steepest_descent(pose)
    ProtoSyn.append(pose, "ilsrr.pdb")
    ProtoSyn.write(pose, "ilsrr_best.pdb")
    outer_ILSRR(pose)
    # inner_monte_carlo(pose)
    println("Seed: $seed")
end