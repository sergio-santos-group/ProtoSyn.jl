push!(LOAD_PATH, "../src")

using ProtoSyn
using ProtoSyn.Peptides

println("ProtoSyn loaded successfully.")

# Example 1.
# -> Create a new segment from an aminoacid string
# There are 2 steps for this task. First, a library of the individual residue
# topologies needs to be loaded. This residues will then be assembled in the
# requested order using the 'build' function. This assembly is performed in
# internal coordinates, same as most methods in ProtoSyn. In order to correctly
# visualize the results, cartesian coordinates need to be generated using the
# 'sync!' function.
#
# Bonus: the @pymol macro allows us to visualize in read-time the results of
# the script. In order to enable the XML-RPC server, launch pymol with the -R
# flag previously. Using this macro without on open server on PyMOL will cause 
# an error.
residue_library = Peptides.loaddb()
pose = Peptides.build("AAA", residue_library)
@pymol sync!(pose)
println("Done buiilding example peptide.")

# QUESTION
# Should the user know when he needs to sync or not? Maybe this task should be
# automatically performed when:
# 1. Building the peptide?
# 2. Visualization functions (@pymol macro or print_to_PDB)?

# Example 2.
# -> Append a single residue to the end of the peptide
# There are 2 steps for this task. First, the required residue needs to be
# loaded as a new fragment from the library, who can then be appended to the
# existing peptide. A fragment is a Pose whose graph is a single segment,
# instead of a full topology. This segment and the original pose segment can
# then be merged. Again, this method is performed in internal coordinates and
# need to be converted to cartesian coordinates using 'sync!' function.
#
# Bonus: PeptideRxToolbelt stands for "peptide reaction toolbelt" and it is a
# set of 3 functions: one to join two peptides, one to split peptides and a
# final function to define the root of the peptide. These are all the functions
# necessary to manipulate peptides and are specific for this type of molecules.
# For example, polycarbohydrates should have different joining, splitting and
# root lookup settings and rules.
# single_residue = Peptides.fragment("P", residue_library)
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