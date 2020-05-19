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
pose = Peptides.build("AAGAAS", residue_library)
@pymol sync!(pose)

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
single_residue = Peptides.fragment("P", residue_library)
@pymol append!(pose, single_residue, 1id, PeptideRxToolbelt)
@pymol @pymol sync!(pose)
