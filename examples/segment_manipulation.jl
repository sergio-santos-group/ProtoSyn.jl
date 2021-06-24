using ProtoSyn

# Load a peptide
pose      = ProtoSyn.Peptides.load("2a3d.pdb")

# Create a peptide from a sequence
res_lib   = ProtoSyn.Peptides.grammar()
appendage = ProtoSyn.Peptides.build(residue_library, seq"ALHYGME",
    ss = ProtoSyn.Peptides.SecondaryStructure[:helix])

# Append both peptides
ProtoSyn.Peptides.append_fragment!(pose, pose.graph[1][end], res_lib, appendage)

# Remove loop
for residue in reverse(rid"29:35"(pose, gather = true))
    ProtoSyn.pop_residue!(pose, residue)
end

# Extrar loop from loaded pose
loop_source_pose = ProtoSyn.Peptides.load("1ctf.pdb")
loop             = fragment(loop_source_pose, rid"40:48")

# Append new loop
ProtoSyn.Peptides.append_fragment!(pose, pose.graph[1][28], res_lib, loop)

# Close the loop with a Monte Carlo simulation
evaluator  = ProtoSyn.Common.default_energy_function()
mutator    = ProtoSyn.Mutators.DihedralMutator(randn, 0.1, 1.0, rid"29:37")
callback   = ProtoSyn.Common.default_energy_step_frame_callback(10, "out.pdb")
thermostat = ProtoSyn.Drivers.get_linear_quench(1.0, 100)

sim = ProtoSyn.Drivers.MonteCarlo(evaluator, mutator, callback, 100, thermostat)
sim(pose)






# ---

















push!(LOAD_PATH, "../src")

using ProtoSyn
using ProtoSyn.Calculators
using ProtoSyn.Peptides
using ProtoSyn
using ProtoSyn.Units
using Printf

# Test
T = Float64
res_lib  = ProtoSyn.Peptides.grammar(T)
sequence   = seq"MRS"
small_pose = ProtoSyn.Peptides.build(res_lib, sequence)
ProtoSyn.write(small_pose, "teste.pdb")

for i in 1:179
    ProtoSyn.ProtoSyn.rotate_dihedral!(small_pose.state, ProtoSyn.Peptides.Dihedral.phi(small_pose.graph[1][2]), deg2rad(2))
    ProtoSyn.append(small_pose, "teste.pdb")
end


pose    = Peptides.load("../2a3d.pdb")
ProtoSyn.write(pose, "teste2.pdb")
ProtoSyn.Peptides.append_fragment!(pose, pose.graph[1][end], res_lib, fragment(small_pose))
ProtoSyn.append(pose, "teste2.pdb")

for residue in reverse(rid"44:50"(pose, gather = true))
    ProtoSyn.Peptides.pop_residue!(pose, residue)
end

ProtoSyn.append(pose, "teste2.pdb")

pose3 = ProtoSyn.Peptides.load("../1c01.pdb")
ProtoSyn.write(pose3, "teste3.pdb")
frag = ProtoSyn.fragment(pose3, rid"35:41")
ProtoSyn.write(Pose(frag), "teste4.pdb")

ProtoSyn.Peptides.append_fragment!(pose, pose.graph[1][43], res_lib, frag)
ProtoSyn.append(pose, "teste2.pdb")

ProtoSyn.bond(pose.graph[1][50]["C"], pose.graph[1][51]["N"])
ProtoSyn.popparent!(pose.graph[1][51]["N"])
ProtoSyn.popparent!(pose.graph[1][51])
ProtoSyn.setparent!(pose.graph[1][51]["N"], pose.graph[1][50]["C"])
ProtoSyn.setparent!(pose.graph[1][51], pose.graph[1][50])
ProtoSyn.append(pose, "teste2.pdb")

energy_function = ProtoSyn.Common.default_energy_function()
energy_function["Bond_Distance_Restraint"].α = 1.0
cb = ProtoSyn.Common.default_energy_step_frame_callback(1, "../teste1.pdb")
ProtoSyn.write(pose, "../teste1.pdb")

selection        = rid"43:50" & an"C$|N$"r
p_mut            = 1/(length(selection(pose, gather = true)))
dihedral_mutator = ProtoSyn.Mutators.DihedralMutator(randn, p_mut, 0.5, selection)
monte_carlo      = ProtoSyn.Drivers.MonteCarlo(
    energy_function, dihedral_mutator, cb, 100,
    ProtoSyn.Drivers.get_linear_quench(0.1, 100))

monte_carlo(pose)

sd = ProtoSyn.Drivers.SteepestDescent(energy_function, cb, 30000, 0.001, 0.1)
sd(pose)



# ---

ProtoSyn.Peptides.remove_sidechains!(pose, res_lib)
ProtoSyn.append(pose, "teste.pdb")
Peptides.cap!(pose)
ProtoSyn.append(pose, "teste.pdb")


# for i in 1:35
#     ProtoSyn.setdihedral!(pose.state, ProtoSyn.Peptides.Dihedral.psi(pose.graph[1][2]), deg2rad(-45 + i * 10))
#     ProtoSyn.append(pose, "teste.pdb")
# end


pose = Peptides.load("../2a3d.pdb")

ef = Calculators.EnergyFunction()
push!(ef, Calculators.TorchANI.get_default_torchani_model(α = 1.0))



pose2 = Peptides.load("../2a3d.pdb")
# rrbm = ProtoSyn.Mutators.RotationRigidBodyMutator(ProtoSyn.rand_vector_in_sphere, randn, ProtoSyn.center_of_mass, 10.0)
# trbm = ProtoSyn.Mutators.TranslationRigidBodyMutator(ProtoSyn.rand_vector_in_sphere, 10.0)
# rrbm(pose2)
# trbm(pose2)
ProtoSyn.write(pose1, "../pose1.pdb")
ProtoSyn.write(pose2, "../pose2.pdb")
ProtoSyn.setdihedral!(pose2.state, Peptides.Dihedral.phi(pose2.graph[1][22]), deg2rad(76))
sync!(pose2)
ProtoSyn.align!(pose2, pose1, an"CA")
ProtoSyn.append(pose1, "../pose1.pdb")
ProtoSyn.append(pose2, "../pose2.pdb")

ProtoSyn.rmsd(pose2, pose1, an"CA")

using ProtoSyn.Peptides
using ProtoSyn
T        = Float64
res_lib  = Peptides.grammar(T)
sequence = seq"MGSWAEFKQRLAAIKTRLQALGGSEAELAAFEKEIAAFESELQAYKGKGNPEVEALRKEAAAIRDELQAYRHN"
pose     = Peptides.build(res_lib, sequence);
begin
    Peptides.setss!(pose, SecondaryStructure[:helix], rid"1:20");
    Peptides.setss!(pose, SecondaryStructure[:helix], rid"27:44");
    Peptides.setss!(pose, SecondaryStructure[:helix], rid"50:end");
    Peptides.remove_sidechains!(pose)
end

pose = ProtoSyn.load("../2a3d.pdb")
ProtoSyn.write(pose, "teste.pdb")

energy_function = ProtoSyn.Common.default_energy_function()
cm_restraints   = Peptides.Calculators.Restraints.load_contact_map("contact_map_example.txt")
energy_function.components[cm_restraints] = 1.0

# pose     = Peptides.build(res_lib, sequence);

selection          = (rid"45" | rid"60") & an"CA"
crankshaft_mutator = ProtoSyn.Mutators.CrankshaftMutator(
    randn, 1.0, 0.5, selection, an"^C$|^O$"r)

ProtoSyn.write(pose, "monte_carlo.pdb")
crankshaft_mutator(pose)
ProtoSyn.append(pose, "monte_carlo.pdb")

e_fcn    = ProtoSyn.Common.default_energy_function()
e_fcn.components[Calculators.Restraints.bond_distance_restraint] = 0.05
cb       = ProtoSyn.Common.default_energy_step_frame_callback(100, "../teste1.pdb")
sd       = ProtoSyn.Drivers.SteepestDescent(e_fcn, cb, 3000, 0.001, 0.1)
brm      = ProtoSyn.Mutators.BlockRotMutator(ProtoSyn.rand_vector_in_sphere, randn, ProtoSyn.center_of_mass, 1.0, 0.2, [rid"1:20", rid"27:44", rid"50:end"], sd)
ProtoSyn.write(pose, "../teste1.pdb")
@time brm(pose)
# ProtoSyn.append(pose, "../teste1.pdb")

begin
    pose = Peptides.build(res_lib, sequence);
    Peptides.setss!(pose, SecondaryStructure[:helix], rid"1:20");
    Peptides.setss!(pose, SecondaryStructure[:helix], rid"27:44");
    Peptides.setss!(pose, SecondaryStructure[:helix], rid"50:end");
end

# rrbm = ProtoSyn.Mutators.RotationRigidBodyMutator(ProtoSyn.rand_vector_in_sphere, randn, ProtoSyn.center_of_mass, 1.0, rid"27:44")
ProtoSyn.write(pose, "../teste1.pdb")
for i in 1:100
    brm(pose)
    ProtoSyn.append(pose, "../teste1.pdb")
end


selection          = (rid"21:26" | rid"45:49") & an"CA" & !rn"PRO"
n                  = count(selection(pose))
p_mut              = 2/(n*(n-1))
crankshaft_mutator = ProtoSyn.Mutators.CrankshaftMutator(
    randn, p_mut, 2.0, selection, an"^C$|^O$"r)


trb_mutator = ProtoSyn.Mutators.TranslationRigidBodyMutator(ProtoSyn.rand_vector_in_sphere, 0.2, rn"CBZ")
rrb_mutator = ProtoSyn.Mutators.RotationRigidBodyMutator(ProtoSyn.rand_vector_in_sphere, randn, 0.2, nothing, rn"CBZ")
compound_mutator = ProtoSyn.Drivers.Compound([trb_mutator, rrb_mutator])


T = Float64
res_lib = grammar(T);

begin
    pose = Peptides.build(res_lib, seq"KKK")
    ProtoSyn.write(pose, "../teste1.pdb")
    Peptides.mutate!(pose, pose.graph[1][1], res_lib, seq"Y")
    ProtoSyn.append(pose, "../teste1.pdb", model = 2)
end

io = open("../selection.pdb", "w"); ProtoSyn.print_selection(io, pose, ((10.0:rn"CBZ") & !(rn"CBZ" | rn"PRO"))(pose)); close(io)

pose = Peptides.load("../mdC.pdb")
ProtoSyn.write(pose, "../teste1.pdb")
design_mutator = Peptides.Mutators.DesignMutator(1.0, res_lib, ((10.0:rn"CBZ") & !(rn"CBZ" | rn"PRO")))
filter!(aa -> aa != 'P', design_mutator.searchable_aminoacids) # Remove proline from searchable aminoacids
for i in 2:6
    design_mutator(pose)
    ProtoSyn.append(pose, "../teste1.pdb")
end

begin
    # pose = Peptides.build(res_lib, seq"MGSWAEFKQRLAAIKTRLQALGGSEAELAAFEKEIAAFESELQAYKGKGNPEVEALRKEAAAIRDELQAYRHN");
    pose = Peptides.build(res_lib, seq"MGSWAEFKQRLAAIKTRLQALGGSEAELAAFEKEIAAFESELQAY");
    Peptides.setss!(pose, SecondaryStructure[:helix], rid"1:20");
    Peptides.setss!(pose, SecondaryStructure[:helix], rid"27:44");
    # Peptides.setss!(pose, SecondaryStructure[:helix], rid"50:end");
    # Peptides.remove_sidechains!(pose)
end

# rl = Peptides.Rotamers.load_dunbrack(T, ProtoSyn.resource_dir * "/Peptides/dunbrack_rotamers.lib")

energy_function = ProtoSyn.Common.get_default_energy_function()
energy_function(pose)


selection = (rid"21:26" | rid"45:49") & an"C$|N$"r
p_mut = 1/(length(selection(pose, gather = true)))
dihedral_mutator1 = ProtoSyn.Mutators.DihedralMutator(randn, p_mut, 0.5, selection)

selection = an"CA" & !rn"PRO"
selection = (rid"21:26" | rid"45:49") & an"CA" & !rn"PRO"
n = count(selection(pose))
p_mut = 2/(n*(n-1))
crankshaft_mutator = ProtoSyn.Mutators.CrankshaftMutator(randn, p_mut, 2.0, selection)

rotamer_blitz = ProtoSyn.Drivers.RotamerBlitz(energy_function, rl, 3, 1)

compound_driver = ProtoSyn.Drivers.Compound([dihedral_mutator2, rotamer_blitz])

function callback_function(pose::Pose, driver_state::ProtoSyn.Drivers.DriverState)
    @printf("STEP %-5d E= %-10.4f AR= %-5.2f%% T= -%7.5f\n", driver_state.step, pose.state.e[:Total], (driver_state.acceptance_count/driver_state.step)*100, driver_state.temperature)
    io = open("../teste1.pdb", "a"); ProtoSyn.write(io, pose); close(io);
end

mc_callback1  = ProtoSyn.Drivers.Callback(callback_function, 500)
monte_carlo1 = ProtoSyn.Drivers.MonteCarlo(energy_function, dihedral_mutator1, mc_callback1, 100000, ProtoSyn.Drivers.get_linear_quench(0.1, 100000))
# monte_carlo2 = ProtoSyn.Drivers.MonteCarlo(energy_function, compound_driver, mc_callback2, 200, ProtoSyn.Drivers.get_constant_temperature(0.005))

# pose = Peptides.build(res_lib, seq"QQQQQQQQQQQQQQQQQQQQQQ");
io = open("../teste1.pdb", "w"); ProtoSyn.write(io, pose); close(io);
monte_carlo1(pose);
# monte_carlo2(pose)

# ----

rm = Peptides.Mutators.RotamerMutator(rl, 1.0, -1, an"CA")
io = open("../teste1.pdb", "w"); ProtoSyn.write(io, pose); close(io);
rm(pose)
io = open("../teste1.pdb", "a"); ProtoSyn.write(io, pose); close(io);


# Peptides.Rotamers.apply!(pose.state, rl["GLN"][-180.0, -180.0][1], pose.graph[1][2])

sync!(pose)
pose.state[1].t += [1.0, 3.0, 4.0]

p_mut = 1/(ProtoSyn.count_residues(pose.graph) * 3)
dihedral_mutator = ProtoSyn.Mutators.DihedralMutator(randn, p_mut, 0.2, an"C$|N$|CA$"r)

energy_function = ProtoSyn.Common.get_default_energy_function()
energy_function(pose, update_forces = true)

io = open("../teste1.pdb", "w"); ProtoSyn.write(io, pose); close(io);

function callback(pose::Pose, driver_state::ProtoSyn.Drivers.DriverState)
    println(driver_state.step)
    io = open("../teste1.pdb", "a"); ProtoSyn.write(io, pose); close(io);
end

mc = ProtoSyn.Drivers.MonteCarlo(energy_function, dihedral_mutator, 100, 0.1)

mc(callback, pose)

using Printf

begin
    pose = Peptides.build(res_lib, seq"MGSWAEFKQRLAAIKTRLQALGGSEAELAAFEKEIAAFESELQAY");
    Peptides.setss!(pose, SecondaryStructure[:helix], rid"1:20");
    Peptides.setss!(pose, SecondaryStructure[:helix], rid"27:44");

    energy_function(pose)
    best = copy(pose)

    function callback(pose::Pose, driver_state::ProtoSyn.Drivers.DriverState)
        global best
        if pose.state.e[:Total] < best.state.e[:Total]
            best = copy(pose)
        end
        if driver_state.step % 50 == 0
            @printf("%5d | Energy: %10.5f | Best Energy: %10.5f | Acceptance Ratio: %10.5f\n", driver_state.step, pose.state.e[:Total], best.state.e[:Total], driver_state.acceptance_count / driver_state.step)
            io = open("../teste1.pdb", "a"); ProtoSyn.write(io, pose); close(io);
        end
    end

    p_mut = 1/(ProtoSyn.count_residues(pose.graph) * 3)
    dihedral_mutator = ProtoSyn.Mutators.DihedralMutator(randn, 1/21, 0.5, an"C$|N$|CA$"r & rid"20:27")
    
    io = open("../teste1.pdb", "w"); ProtoSyn.write(io, pose); close(io);
    for temp in [0.01, 0.008, 0.005, 0.003, 0.001, 0.0005, 0.0001, 0.00001]
        println(" Temperature: $temp")
        mc = ProtoSyn.Drivers.MonteCarlo(energy_function, dihedral_mutator, 2000, temp)
        mc(callback, pose)
        pose = copy(best)
    end

    sd = ProtoSyn.Drivers.SteepestDescent(energy_function, 100, 0.001)
    sd(callback, pose)
end


io = open("../teste1.pdb", "w"); ProtoSyn.write(io, pose); close(io);
mask = an"C$|N$|CA$"r(pose)
for i in 1:100
    dihedral_mutator(pose, mask)
    sync!(pose)
    io = open("../teste1.pdb", "a"); ProtoSyn.write(io, pose); close(io);
end

f = Calculators.get_energy_function()
f(original)

# results = Dict{String, Vector{Float64}}()
for aminoacid in ["A", "R", "N", "D", "C", "E", "Q", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]
    pose = Peptides.build(res_lib, [aminoacid])

    results = "$aminoacid"
    results *= @sprintf(" %10.5f", ProtoSyn.Calculators.TorchANI.calc_torchani_ensemble(pose))
    for i in 0:7
        results *= @sprintf(" %10.5f", ProtoSyn.Calculators.TorchANI.calc_torchani_model(pose, i))
    end

    println(results)
end


println("ProtoSyn loaded successfully.")

# Example 1.
# -> Create a new segment from an aminoacid string
T = Float64
n_samples = 10
res_lib = grammar(T);
pose = Peptides.build(res_lib, seq"GG");

verlet_list = ProtoSyn.Calculators.VerletList(pose.state.size);
verlet_list.cutoff = 3.0;
@time ProtoSyn.Calculators.update_simd!(verlet_list, pose);
ProtoSyn.Calculators.serial(pose.state, verlet_list)
# pose = Peptides.load("../2a3d.pdb");

ProtoSyn.Calculators.cuda(pose.state);
@time begin
    for i in 1:n_samples
        ProtoSyn.Calculators.cuda(pose.state);
    end
end

verlet_list = ProtoSyn.Calculators.VerletList(pose.state.size);
verlet_list.cutoff = 12.0; # In Angstrom (Å)
@time ProtoSyn.Calculators.update_simd!(verlet_list, pose.state);
@time ProtoSyn.Calculators.update_simd!(verlet_list, pose.state);

ProtoSyn.Calculators.simd(pose.state, verlet_list);
@time begin
    for i in 1:n_samples
        ProtoSyn.Calculators.simd(pose.state, verlet_list);
    end
end

# @time ProtoSyn.Calculators.simd(pose.state, verlet_list, Vec{8, Float64});
# @time ProtoSyn.Calculators.opa(pose.state, verlet_list);

verlet_list = ProtoSyn.Calculators.VerletList(pose.state.size);
verlet_list.cutoff = 12.0;
@time ProtoSyn.Calculators.update_serial!(verlet_list, pose.state);
@time ProtoSyn.Calculators.update_serial!(verlet_list, pose.state);

ProtoSyn.Calculators.serial(pose.state, verlet_list);
@time begin
    for i in 1:n_samples
        ProtoSyn.Calculators.serial(pose.state, verlet_list);
    end
end

# ProtoSyn.Calculators.serial(pose.state)

# julia> ProtoSyn.Calculators.naive(pose.state, Inf)
# 7×7 Array{Float64,2}:
#  0.0  1.00939  1.44898  2.08368  2.08368  2.44014  2.66707
#  0.0  0.0      2.12288  2.83006  2.83006  2.51674  2.24842
#  0.0  0.0      0.0      1.08987  1.08987  1.52233  2.39287
#  0.0  0.0      0.0      0.0      1.78     2.14098  3.08794
#  0.0  0.0      0.0      0.0      0.0      2.14098  3.08794
#  0.0  0.0      0.0      0.0      0.0      0.0      1.22933
#  0.0  0.0      0.0      0.0      0.0      0.0      0.0


# Peptides.append_residues!(pose, pose.graph[1][end], res_lib, seq"GGG");
# unbond!(pose, pose.graph[1][2]["C"], pose.graph[1][3]["N"]);
# sync!(pose);
# println("Outside angle: $(pose.state[pose.graph[1][3]["N"]].ϕ) $(ProtoSyn.dihedral(pose.state[pose.graph[1][3]["N"]], pose.state[ProtoSyn.root(pose.graph)], pose.state[ProtoSyn.root(pose.graph).parent], pose.state[ProtoSyn.root(pose.graph).parent.parent]))")
# Peptides.insert_residues!(pose, pose.graph[1][1], res_lib, seq"LLL");
# Peptides.mutate!(pose, pose.graph[1][2], res_lib, seq"L");
# setss!(pose, SecondaryStructure[:helix], rid"2:end");

# pose = Peptides.build(res_lib, seq"AAAAAAAAAA");
# ProtoSyn.setdihedral!(pose.state, pose.graph[1][3]["CA"], 0°); # Bug: is measuring i-2 and setting i-1
# io = open("../teste1.pdb", "w"); ProtoSyn.write(io, pose); close(io);
# pose = Peptides.load("../teste1.pdb");
# ProtoSyn.setdihedral!(pose.state, pose.graph[1][3]["CA"], 180°); # Bug: is measuring i-2 and setting i-1
# io = open("../teste2.pdb", "w"); ProtoSyn.write(io, pose); close(io);

# Peptides.setdihedral!(pose.state, pose.graph[1][23], Peptides.Dihedral.phi, -60°)
# Peptides.setdihedral!(pose.state, pose.graph[1][23], Peptides.Dihedral.psi, 150°)
# Peptides.setdihedral!(pose.state, pose.graph[1][2], Peptides.Dihedral.omega, 180°)
# ProtoSyn.setdihedral!(pose.state, pose.graph[1][3]["C"], 0°);
# Peptides.setss!(pose, SecondaryStructure[:helix], rid"2:end");

# pose = Peptides.load("../2a3d.pdb");

using PyCall
torch    = pyimport("torch")
torchani = pyimport("torchani")

if torch.cuda.is_available()
    device = torch.device("cuda")
else
    device = torch.device("cpu")
end

model = torchani.models.ANI2x(periodic_table_index = true).to(device)

function calc_ani_energy(pose::Pose, model)
    c = get_cartesian_matrix(pose)
    coordinates = torch.tensor([c], requires_grad = true, device = device).float()

    s = get_ani_species(pose)
    species = torch.tensor([s], device = device)

    return float(model((species, coordinates)))
end

function model2(t::Tuple)

    global model

    m1 = model.species_converter(t)
    m2 = model.aev_computer(m1)
    m3 = get(model.neural_networks, 0)(m2)
    # println(m3)
    return m3[2].item()
    # return model.energy_shifter(m3)
end

# ---

function uncap(pose) # Pyrosetta version

    ris = pyrosetta.rosetta.core.select.residue_selector.ResidueIndexSelector
    cap_indexes = ""
    for residue in pose.residues
        if residue.is_terminus()
            cap_indexes *= "$(residue.seqpos()),"
        end
    end
    println(cap_indexes[1:(end-1)])
    caps = ris(cap_indexes[1:(end-1)])

    uncapper = pyrosetta.rosetta.protocols.simple_moves.ModifyVariantTypeMover()
    uncapper.set_additional_type_to_remove("LOWER_TERMINUS_VARIANT")
    uncapper.set_additional_type_to_remove("UPPER_TERMINUS_VARIANT")
    uncapper.set_residue_selector(caps)
    uncapper.set_update_polymer_bond_dependent_atoms(true)

    uncapper.apply(pose)
end


pyrosetta = pyimport("pyrosetta")
pyrosetta.init()

pyrosetta_score_fn = pyrosetta.get_fa_scorefxn()

function calc_pyrosetta_energy(pose::Pose)

    global pyrosetta_score_fn
    global pyrosetta

    pyrosetta_pose = pyrosetta.pose_from_sequence(Peptides.get_sequence(pose))
    # Set coords TO DO
    return pyrosetta_score_fn(pyrosetta_pose)
end

# ---

pose1 = Peptides.build(res_lib, seq"KKKKKKKKKKKKKKKKK");
pose = Peptides.build(res_lib, seq"NNNNNNNNNNNNNNNNNN");
pose = Peptides.load("../2a3d.pdb");

n_steps = 1000
t = 0.005
print_every = 50
step_size = 0.15

e = calc_ani_energy(pose, model2)
saved_state = copy(pose.state)
accepted = 0

io = open("../teste.pdb", "w"); ProtoSyn.write(io, pose); close(io);

for step in 1:n_steps
    global e
    global saved_state
    global accepted
    global print_every
    global step_size
    
    r = rand(pose.graph.items[1].items)
    phi_psi = rand([Peptides.Dihedral.phi, Peptides.Dihedral.psi])
    val = randn() * step_size # in radians
    Peptides.rotate_dihedral!(pose.state, r, phi_psi, val)
    ne = calc_ani_energy(pose, model2)
    n = rand()

    if (ne < e) || (n < exp((e - ne)/ t))
        e = ne
        saved_state = copy(pose.state)
        accepted += 1
    else
        pose.state = copy(saved_state)
    end

    if step % print_every == 0
        @printf("Step: %4d | Energy: %10.3f | Acceptance Rate: %s\n", step, e, accepted/step) 
        io = open("../teste.pdb", "a"); ProtoSyn.write(io, pose); close(io);
    end
end

calc_ani_energy(pose, model2)



# residues = rid"1:end"(pose, gather = true)
# for r in residues
#     Peptides.setdihedral!(pose.state, r, Peptides.Dihedral.phi, -60°)
#     Peptides.setdihedral!(pose.state, r, Peptides.Dihedral.psi,  -45°)
#     Peptides.setdihedral!(pose.state, r, Peptides.Dihedral.omega, 180°)
# end
Peptides.setss!(pose, SecondaryStructure[:helix], rid"2:end");
io = open("../teste111.pdb", "w"); ProtoSyn.write(io, pose); close(io);

pose = Peptides.load("../2a3d.pdb");
Peptides.uncap!(pose, pose.graph[1][1]);
io = open("../teste2.pdb", "w"); ProtoSyn.write(io, pose); close(io);
Peptides.insert_residues!(pose, pose.graph[1][1], res_lib, seq"LLL");
io = open("../teste3.pdb", "w"); ProtoSyn.write(io, pose); close(io);
# Peptides.unbond!(pose, pose.graph[1][4], pose.graph[1][5]);
# Peptides.setss!(pose, SecondaryStructure[:linear], pose.graph[1].items[1:4]);
# Peptides.append_residues!(pose, pose.graph[1][4], res_lib, seq"LLL");

# ProtoSyn.setdihedral!(pose.state, pose.graph[1][3]["CA"], 0°);
# Peptides.setdihedral!(pose.state, pose.graph[1][3], Peptides.Dihedral.omega, 0°);
# ProtoSyn.setdihedral!(pose.state, pose.graph[1][7]["CA"], 0°);
# pop!(pose, pose.graph[1][9])





begin
    @printf("%-40s > %8s | %13s | %8s | %8s\n", "Atom", "b", "θ", "ϕ", "Δϕ")
    for atom in ProtoSyn.eachatom(pose.graph)
        @printf("%-40s > %8.3f | %13.8f | %8.3f | %8.3f\n", atom, pose.state[atom].b, rad2deg(pose.state[atom].θ), rad2deg(pose.state[atom].ϕ), rad2deg(pose.state[atom].Δϕ))
    end
end

# s = pose.state
# a = pose.graph[1]
# ProtoSyn.request_i2c!(pose.state)
# sync!(pose)

# println("PHI 22: $(rad2deg(ProtoSyn.dihedral(s[a[21]["C"]], s[a[22]["N"]], s[a[22]["CA"]], s[a[22]["C"]])))")
# for atom in [Peptides.Dihedral.phi, Peptides.Dihedral.psi, Peptides.Dihedral.omega]
#     println("Dihedral 22 (INTERNAL): $(rad2deg(Peptides.getdihedral(pose, a[22], atom))) (atom $(atom)-$(a[22][atom.atom].ascendents))")
# end
# pose = Peptides.build(res_lib, seq"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA");
io = open("../teste.pdb", "w"); ProtoSyn.write(io, pose); close(io);

res_lib = grammar();
pose = Peptides.build(res_lib, seq"AAAA");
Peptides.setdihedral!(pose, pose.graph[1][2], Peptides.Dihedral.phi, deg2rad(90));

# unbond!(pose, pose.graph[1][2]["C"], pose.graph[1][3]["N"]);
# ProtoSyn.bond(pose, pose.graph[1][2]["C"], pose.graph[1][3]["N"], res_lib);
# pop!(pose, pose.graph[1][1])
# Peptides.insert_residues!(pose, pose.graph[1][1], res_lib, seq"KKK");
# Peptides.insert_residues!(pose, pose.graph[1][3], res_lib, seq"Q");
# Peptides.append_residues!(pose, pose.graph[1][3], res_lib, seq"LLL");

# Peptides.insert_residues!(pose, pose.graph[1][3], res_lib, seq"Q");
# Peptides.insert_residues!(pose, pose.graph[1][8], res_lib, seq"Q");
# Peptides.insert_residues!(pose, pose.graph[1][10], res_lib, seq"Q");
# setdihedral!(pose.state, pose.graph[1][2][Peptides.DihedralType.phi], deg2rad(-90));

# Peptides.append_residues!(pose, pose.graph[1][end], res_lib, seq"GGG");
Peptides.insert_residues!(pose, pose.graph[1][6], res_lib, seq"G");





# When prepending to a severed segment, since the new residues added will be
# connected to the origin, there's no way for ProtoSyn to know where to place
# the new residues, and therefore they will show overlapping the origin
# sync!(pose)
# println(pose.state[pose.graph[1][6]["N"]].b)
# first_res = pose.graph[1][1];
# Peptides.prepend_residues!(pose, first_res, res_lib, seq"GG", ss = SecondaryStructure[:linear], op = "α");
# println(pose.state[pose.graph[1][8]["N"]].b)

# Note that, since residues 6 and 7 have been unbonded, residue 7's origin is
# now the global topology origin. Therefore, if setting the phi angle on residue
# 7, this will be performed in relation to the origin position, leading to major
# changes in position. Therefore, when setting secondary structures, the
# origin.children residues should not changed.
setss!(pose, SecondaryStructure[:helix]);


first_res = pose.graph[1][1];
Peptides.prepend_residues!(pose, first_res, res_lib, seq"GGG", ss = SecondaryStructure[:helix], op = "α");

last_res = pose.graph[1][end];
Peptides.append_residues!(pose, last_res, res_lib, seq"KKK", ss = SecondaryStructure[:parallel_sheet], op = "α");



setss!(pose, SecondaryStructure[:linear])

# Example 2.
# -> Append a single residue to the end of the peptide, joining it to the existing segment
Peptides.append_residues!(pose, last_res, res_lib, seq"AAAAAAAAAAAAAAA", ss = SecondaryStructure[:linear], op = "α")


# ProtoSyn.join_all_segments(pose)
# setss!(pose, SecondaryStructure[:linear], MaxSerialSelection{Segment}(:index))
# sync!(pose)

# @pymol append_fragment_as_new_segment!(pose, single_residue, 1id, PeptideRxToolbelt)
# @pymol sync!(pose)

# Optinally, one can use a wrapping function to append a residue to the end
# of an existing pose.
# @pymol append_fragment_as_new_segment!(pose, "A")

# Example 3.
# -> Append multiple residues to the end of the peptide
# The same functions used in the previous example can be employed to add more
# than 1 residue simultaneously to the end of the petidic chain.
# multiple_residues = Peptides.fragment("PAPA", residue_library)
# @pymol append_fragment_as_new_segment!(pose, multiple_residues, 1id, PeptideRxToolbelt)
# @pymol sync!(pose)

@pymol append_fragment_as_new_segment!(pose, "PAPA")

# QUESTION
# When adding mutiple residues, the resulting conformation is not linear, it has
# some deviation

# QUESTION
# Selections always return atom arrays?

# Example 4.
# -> Remove the last residue on the peptidic segment
# The 
pop!(pose, pose.graph[1][end])
@pymol sync!(pose)

# QUESTION
# When doing @pymol pop!(pose, pose.graph[1][end]) it gives and error:
# Attempt to access 10-element Array{AtomState{Float64},1} at index [69] ?

# Example 5.
# -> Remove a residue in the center of the peptidic segment
# The same pp
pop!(pose, pose.graph[1, 2])
@pymol sync!(pose)