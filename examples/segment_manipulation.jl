push!(LOAD_PATH, "../src")

using ProtoSyn
using ProtoSyn.Peptides
using ProtoSyn.Builder

println("ProtoSyn loaded successfully.")

# Example 1.
# -> Create a new segment from an aminoacid string
res_lib = grammar();
pose = Builder.build(res_lib, seq"AAAAQR");
io = open("../teste.pdb", "w"); ProtoSyn.write(io, pose); close(io)

setss!(pose, SecondaryStructure[:linear])
sync!(pose)

# Example 2.
# -> Append a single residue to the end of the peptide, joining it to the existing segment
last_res = pose.graph[1][end]
Peptides.append_residues(pose, last_res, res_lib, seq"AAAAAAAAAAAAAAA", ss = :helix, op = "α")


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