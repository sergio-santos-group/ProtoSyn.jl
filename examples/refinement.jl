using ProtoSyn
using ProtoSyn.Peptides
using ProtoSyn
using Printf

# 1) Load the default residue library
T = Float64
res_lib = grammar(T)

# 2) Load the default Dunbrack rotamer library
rot_lib = Peptides.Rotamers.load_dunbrack(T)

# 3) Load the initial pose (and add sidechains)
pose = Peptides.load("monte_carlo_1.pdb")
ProtoSyn.Peptides.add_sidechains!(pose, res_lib)
ProtoSyn.write(pose, "refinement.pdb")

# 5) Load the custom energy functions
energy_function_1 = ProtoSyn.Calculators.EnergyFunction(Dict(
    ProtoSyn.Calculators.TorchANI.torchani_model => T(1.0),
    ProtoSyn.Peptides.Calculators.Restraints.sidechain_clash_restraint => T(0.1),
    ProtoSyn.Calculators.Restraints.bond_distance_restraint => T(10000.0)
))

energy_function_2 = ProtoSyn.Calculators.EnergyFunction(Dict(
    ProtoSyn.Calculators.TorchANI.torchani_model => T(1.0)
))

# 6) Define the RotamerBlitz driver (Here we will consider the 3 most likely
# rotamers, with 1 pass). Note that only the TorchANI component is being
# employed (with sidechain clash it wouldbe too slow), and that Proline residues
# won't be sampled.
rotamer_blitz = Peptides.Drivers.RotamerBlitz(energy_function_2, rot_lib, 3, 1)
rotamer_blitz.selection = !rn"PRO"

# 7) Define the other mutators of the first stage (Dihedral and Crankshaft)
# In this example, no blocks are being kept artificially in the secondary
# structure, so small bends and kinks are possible. However, the pmut and step
# size values have been lowered. 
selection          = an"C$|N$"r
p_mut              = 1/(count(selection(pose)))
dihedral_mutator   = ProtoSyn.Mutators.DihedralMutator(
    randn, p_mut, 0.05, selection)

selection          = an"CA" & !rn"PRO"
n                  = count(selection(pose))
p_mut              = 0.1/(n*(n-1))
crankshaft_mutator = ProtoSyn.Mutators.CrankshaftMutator(
    randn, p_mut, 0.005, selection, !(an"^CA$|^N$|^C$|^H$|^O$"r))

# 7) Define the dihedral + crankshaft compound driver
compound_driver_1 = ProtoSyn.Drivers.CompoundDriver(
    [dihedral_mutator, crankshaft_mutator, rotamer_blitz])

# 8) Define the stage 1 callback
function get_callback_function(stage::Int)
    return function callback_function(pose::Pose, driver_state::ProtoSyn.Drivers.DriverState)
        s = @sprintf(" STAGE_%d STEP %-6d | E(total)= %-10.4f | AR= %-5.1f%%  | T= %-7.4f | E(ani)= %-10.4f | E(cls)= %-10.4f\n", stage, driver_state.step, pose.state.e[:Total], (driver_state.acceptance_count/driver_state.step)*100, driver_state.temperature, pose.state.e[:TorchANI_ML_Model], pose.state.e[:Clash_Sidechan_Restraint])
        print(s)
        open("refinement.dat", "a") do io
            write(io, s)
        end
        driver_state.step == 0 && return
        ProtoSyn.append(pose, "refinement.pdb")
    end
end
callback_1 = ProtoSyn.Drivers.Callback(get_callback_function(1), 1)

# 9) Define the Stage 1 Monte Carlo Driver
T_i_1         = 0.01
n_steps_1     = 300
monte_carlo_1 = ProtoSyn.Drivers.MonteCarlo(
    energy_function_1,
    compound_driver_1,
    callback_1,
    n_steps_1,
    ProtoSyn.Drivers.get_linear_quench(T_i_1, n_steps_1))

# 10) Run the Stage 1 Monte Carlo Driver
energy_function_1(pose)
s = @sprintf(" STAGE_1 STEP %-6d | E(total)= %-10.4f | AR= %-5.1f%%  | T= %-7.4f | E(ani)= %-10.4f | E(cls)= %-10.4f\n", 0, pose.state.e[:Total], 0.0, 0.0, pose.state.e[:TorchANI_ML_Model], pose.state.e[:Clash_Sidechan_Restraint])
print(s)
open("refinement.dat", "w") do io
    write(io, s)
end
ProtoSyn.ProtoSyn.write(pose, "refinement.pdb")
monte_carlo_1(pose)

# ---

# 11) Load the Stage 2 Mutators (BlockRot + BlockTrans)
cb = ProtoSyn.Common.default_energy_step_callback(100)
sd = ProtoSyn.Drivers.SteepestDescent(energy_function_2, cb, 1000, 0.001, 0.1)
brm = ProtoSyn.Mutators.BlockRotMutator(ProtoSyn.rand_vector_in_sphere, randn, ProtoSyn.center_of_mass, 0.3, 0.3, [rid"1:20", rid"27:44", rid"50:end"], nothing)
btm = ProtoSyn.Mutators.BlockTransMutator(ProtoSyn.rand_vector_in_sphere, 0.3, 0.3, [rid"1:20", rid"27:44", rid"50:end"], sd)
blockrot_driver = ProtoSyn.Drivers.CompoundDriver([brm, btm])

compound_driver_2 = ProtoSyn.Drivers.CompoundDriver(
    [blockrot_driver, rotamer_blitz])

# 12) Define the Stage 1 Monte Carlo Driver
T_i_2         = 0.05
n_steps_2     = 200
callback_2    = ProtoSyn.Drivers.Callback(get_callback_function(2), 1)
monte_carlo_2 = ProtoSyn.Drivers.MonteCarlo(
    energy_function_1,
    compound_driver_2,
    callback_2,
    n_steps_2,
    ProtoSyn.Drivers.get_linear_quench(T_i_2, n_steps_2))

# 13) Run the Stage 1 Monte Carlo Driver
monte_carlo_2(pose)
# ---

# 11) Run the steepest descent
cb = ProtoSyn.Common.default_energy_step_frame_callback(5, "refinement.pdb")
sd = ProtoSyn.Drivers.SteepestDescent(energy_function_rb, cb, 1000, 0.001, 0.01)

sd(pose)

base_truth = ProtoSyn.Peptides.load("../2a3d.pdb")
energy_function_rb(base_truth)


# ---

e = energy_function(pose)
println(" (before) Energy: $e")
println("Starting rotamer blitz ...")
rotamer_blitz(pose)
e = energy_function(pose)
println(" (after) Energy: $e")
ProtoSyn.append(pose, "refinement.pdb")

# ---


function progress_bar(bars::Vector{Tuple{Int, Int}})
    print("\u1b[2F")
    for bar in bars
        current, total = bar[1], bar[2]
        @assert current <= total "Can't write a progress bar over 100% full."
        c = floor(Int, (current/total)*100)
        print("|"*"#"^c*"-"^(100-c)*"|\n")
    end
    print("\u1b[2K")
    println()
end

pose = ProtoSyn.Peptides.load("monte_carlo_1.pdb")
ProtoSyn.Peptides.add_sidechains!(pose, res_lib)
ProtoSyn.write(pose, "refinement.pdb")

energy_function = ProtoSyn.Common.default_energy_function()
energy_function.components[ProtoSyn.Peptides.Calculators.Restraints.sidechain_clash_restraint] = 1.0
energy_function.components[ProtoSyn.Calculators.Restraints.bond_distance_restraint] = 100.0

sequence = seq"MGM"
pose     = ProtoSyn.Peptides.build(res_lib, sequence)
ProtoSyn.setdihedral!(pose.state, ProtoSyn.Peptides.Dihedral.chi1(pose.graph[1][3]), rad2deg(-10))
ProtoSyn.write(pose, "refinement_t.pdb")

cb = ProtoSyn.Common.default_energy_step_frame_callback(5, "refinement.pdb")
sd = ProtoSyn.Drivers.SteepestDescent(energy_function, cb, 1000, 0.001, 0.01)

sd(pose)
sync!(pose)

for rid in 1:1 #length(eachresidue(pose.graph))
    sele = SerialSelection{Residue}(rid, :id)
    ProtoSyn.Peptides.add_sidechains!(pose, res_lib)

    ProtoSyn.write(pose, "refinement.pdb")

    energy_function(pose, update_forces = true)
    sd(pose)
end

# ---
exit(1)

function print_forces(pose)
    open("teste.forces", "w") do file_out
        for (i, atom) in enumerate(eachatom(pose.graph))
            !any(k -> k != 0, pose.state.f[:, i]) && continue
            x  = pose.state[atom].t[1]
            y  = pose.state[atom].t[2]
            z  = pose.state[atom].t[3]
            fx = x + pose.state.f[1, i]
            fy = y + pose.state.f[2, i]
            fz = z + pose.state.f[3, i]
            s  = @sprintf("%5d %12.3f %12.3f %12.3f %12.3f %12.3f %12.3f\n", atom.id, x, y, z, fx, fy, fz)
            write(file_out, s)
        end
    end
end

α_sol   = 0.01
α_cmap  = 0.001
T_i     = 0.1
n_steps = 100

energy_function = ProtoSyn.Common.default_energy_function()
cm_restraints   = Peptides.Calculators.Restraints.load_contact_map("contact_map_example.txt")
energy_function.components[ProtoSyn.Peptides.Calculators.Caterpillar.solvation_energy] = α_sol
energy_function.components[cm_restraints] = α_cmap
energy_function(pose)

sd_energy_function = copy(energy_function)
sd_energy_function.components[cm_restraints] = 0.0
cb = ProtoSyn.Common.default_energy_step_callback(500)
sd = ProtoSyn.Drivers.SteepestDescent(sd_energy_function, cb, 1000, 0.001, 0.1)

brm = ProtoSyn.Mutators.BlockRotMutator(ProtoSyn.rand_vector_in_sphere, randn, ProtoSyn.center_of_mass, 0.3, 0.3, [rid"1:20", rid"27:44", rid"50:end"], nothing)
btm = ProtoSyn.Mutators.BlockTransMutator(ProtoSyn.rand_vector_in_sphere, 0.3, 0.3, [rid"1:20", rid"27:44", rid"50:end"], sd)

compound_driver = ProtoSyn.Drivers.CompoundDriver([brm, btm])

function callback_function(pose::Pose, driver_state::ProtoSyn.Drivers.DriverState)
    s = @sprintf(" INNER STEP %-6d | E(total)= %-10.4f | AR= %-5.1f%%  | T= %-7.4f | E(ani)= %-10.4f | E(sol)= %-10.4f | E(cmp)= %-10.4f | E(cls)= %-10.4f\n", driver_state.step, pose.state.e[:Total], (driver_state.acceptance_count/driver_state.step)*100, driver_state.temperature, pose.state.e[:TorchANI_ML_Model], pose.state.e[:Caterpillar_Solvation], pose.state.e[:Contact_Map_Restraint], pose.state.e[:Clash_Restraint])
    print(s)
    open("monte_carlo.dat", "a") do io
        write(io, s)
    end
    driver_state.step == 0 && return
    ProtoSyn.append(pose, "monte_carlo.pdb")
end

callback    = ProtoSyn.Drivers.Callback(callback_function, 1)
monte_carlo = ProtoSyn.Drivers.MonteCarlo(
    energy_function,
    compound_driver,
    callback,
    n_steps,
    ProtoSyn.Drivers.get_linear_quench(T_i, n_steps))

ProtoSyn.write(pose, "monte_carlo.pdb")
monte_carlo(pose)