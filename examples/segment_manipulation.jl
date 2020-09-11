push!(LOAD_PATH, "../src")

using ProtoSyn
using ProtoSyn.Peptides
using ProtoSyn.Builder
using ProtoSyn.Units

println("ProtoSyn loaded successfully.")

# Example 1.
# -> Create a new segment from an aminoacid string

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

# Peptides.setdihedral!(pose.state, pose.graph[1][2], Peptides.Dihedral.phi, -60°)
# Peptides.setdihedral!(pose.state, pose.graph[1][2], Peptides.Dihedral.psi, -45°)
# Peptides.setdihedral!(pose.state, pose.graph[1][2], Peptides.Dihedral.omega, 180°)
# Builder.setdihedral!(pose.state, pose.graph[1][3]["C"], 0°);
# Peptides.setss!(pose, SecondaryStructure[:helix], rid"2:end");

res_lib = grammar();
pose = Peptides.build(res_lib, seq"AAAAAAAAAAAA");
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