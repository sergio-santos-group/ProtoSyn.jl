push!(LOAD_PATH, "../src")

using ProtoSyn
using ProtoSyn.Calculators
using ProtoSyn.Peptides
using ProtoSyn.Builder
using ProtoSyn.Units

res_lib = grammar(Float64);
pose = Peptides.build(res_lib, seq"GGGGGGGGGGGGGGGGGGGG");

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
# unbond(pose, pose.graph[1][2]["C"], pose.graph[1][3]["N"]);
# sync!(pose);
# println("Outside angle: $(pose.state[pose.graph[1][3]["N"]].ϕ) $(ProtoSyn.dihedral(pose.state[pose.graph[1][3]["N"]], pose.state[ProtoSyn.origin(pose.graph)], pose.state[ProtoSyn.origin(pose.graph).parent], pose.state[ProtoSyn.origin(pose.graph).parent.parent]))")
# Peptides.insert_residues!(pose, pose.graph[1][1], res_lib, seq"LLL");
# Peptides.mutate!(pose, pose.graph[1][2], res_lib, seq"L");
# setss!(pose, SecondaryStructure[:helix], rid"2:end");

# pose = Peptides.build(res_lib, seq"AAAAAAAAAA");
# Builder.setdihedral!(pose.state, pose.graph[1][3]["CA"], 0°); # Bug: is measuring i-2 and setting i-1
# io = open("../teste1.pdb", "w"); ProtoSyn.write(io, pose); close(io);
# pose = Peptides.load("../teste1.pdb");
# Builder.setdihedral!(pose.state, pose.graph[1][3]["CA"], 180°); # Bug: is measuring i-2 and setting i-1
# io = open("../teste2.pdb", "w"); ProtoSyn.write(io, pose); close(io);

# Peptides.setdihedral!(pose.state, pose.graph[1][23], Peptides.Dihedral.phi, -60°)
# Peptides.setdihedral!(pose.state, pose.graph[1][23], Peptides.Dihedral.psi, 150°)
# Peptides.setdihedral!(pose.state, pose.graph[1][2], Peptides.Dihedral.omega, 180°)
# Builder.setdihedral!(pose.state, pose.graph[1][3]["C"], 0°);
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
# Peptides.unbond(pose, pose.graph[1][4], pose.graph[1][5]);
# Peptides.setss!(pose, SecondaryStructure[:linear], pose.graph[1].items[1:4]);
# Peptides.append_residues!(pose, pose.graph[1][4], res_lib, seq"LLL");

# Builder.setdihedral!(pose.state, pose.graph[1][3]["CA"], 0°);
# Peptides.setdihedral!(pose.state, pose.graph[1][3], Peptides.Dihedral.omega, 0°);
# Builder.setdihedral!(pose.state, pose.graph[1][7]["CA"], 0°);
# pop!(pose, pose.graph[1][9])





begin
    @printf("%-40s > %8s | %13s | %8s | %8s\n", "Atom", "b", "θ", "ϕ", "Δϕ")
    for atom in ProtoSyn.eachatom(pose.graph)
        @printf("%-40s > %8.3f | %13.8f | %8.3f | %8.3f\n", atom, pose.state[atom].b, rad2deg(pose.state[atom].θ), rad2deg(pose.state[atom].ϕ), rad2deg(pose.state[atom].Δϕ))
    end
end

# s = pose.state
# a = pose.graph[1]
# ProtoSyn.request_i2c(pose.state)
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

# unbond(pose, pose.graph[1][2]["C"], pose.graph[1][3]["N"]);
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

# @pymol append!(pose, single_residue, 1id, PeptideRxToolbelt)
# @pymol sync!(pose)

# Optinally, one can use a wrapping function to append a residue to the end
# of an existing pose.
# @pymol append!(pose, "A")

# Example 3.
# -> Append multiple residues to the end of the peptide
# The same functions used in the previous example can be employed to add more
# than 1 residue simultaneously to the end of the petidic chain.
# multiple_residues = Peptides.fragment("PAPA", residue_library)
# @pymol append!(pose, multiple_residues, 1id, PeptideRxToolbelt)
# @pymol sync!(pose)

@pymol append!(pose, "PAPA")

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