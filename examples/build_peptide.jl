# In this example we will explore how to build and manipulate a peptidic
# system (appending and removing residues, etc).
# To build and manipulate peptides, we must first load a residue library, also
# know in ProtoSyn as a `grammar`. This grammar, besides the default topologies
# of aminoacids, also contains instructions on how to connect two aminoacids
# (the peptide bond), etc. The Peptides module of ProtoSyn makes available a
# default grammar with the 20 natural aminoacids.
# Note: Most functions used in this example have a counterpart in base ProtoSyn.
# However, when dealing with peptides, the Peptide module functions should be
# used.

using ProtoSyn
using ProtoSyn
using ProtoSyn.Peptides
using ProtoSyn.Units

T = Float64
res_lib = Peptides.grammar(T)

# Using this grammar, we can now build a peptide in the default linear
# conformation
pose = Peptides.build(res_lib, seq"AGGQMKG")

# Note: We can specify a SecondaryStructure when building a new peptide.
pose = Peptides.build(res_lib, seq"AGQMK", Peptides.SecondaryStructure[:helix])

# We can visualize our newly created peptide by printing to PDB format:
ProtoSyn.write(pose, "build_peptide.pdb")

# ------------------------------------------------------------------------------
# (A) Changing the secondary structure
Peptides.setss!(pose, Peptides.SecondaryStructure[:antiparallel_sheet])

# We can append new frames to an existing PDB file using the ProtoSyn.append 
# function:
ProtoSyn.append(pose, "build_peptide.pdb", model = 2)

# ------------------------------------------------------------------------------
# (B) Appending one or more residues to the end of a structure. 
Peptides.append_residues!(pose, pose.graph[1][end], res_lib, seq"GGG");

# We can also append the residues in a given secondary structure:
Peptides.append_residues!(pose, pose.graph[1][end], res_lib, seq"GGG",
    ss = Peptides.SecondaryStructure[:helix])

# ------------------------------------------------------------------------------
# (C) Remove residues from the end of a structure.
Peptides.pop_residue!(pose, pose.graph[1][end])

# ------------------------------------------------------------------------------
# (D) Insert residues in the beggining/middle of the structure
# Note: This function inserts the residue on the same position as the given
# residue, but does not maintain the secondary structure. The phi, psi and omega
# angles can, however, be supplied in optional `ss` parameter. You can also add
# multiple residues in a row.
Peptides.insert_residues!(pose, pose.graph[1][1], res_lib, seq"K")
Peptides.insert_residues!(pose, pose.graph[1][3], res_lib, seq"KK")

# ------------------------------------------------------------------------------
# (E) Set secondary structure of the whole structure or selection
# Note: See the selection.jl example for more info on selections.
Peptides.setss!(pose, Peptides.SecondaryStructure[:helix])
Peptides.setss!(pose, Peptides.SecondaryStructure[:helix], rn"LYS")

# ------------------------------------------------------------------------------
# (F) Mutate an aminoacid
# Note: the `mutate!` function maintains the pre-existent secondary structure.
Peptides.mutate!(pose, pose.graph[1][3], res_lib, seq"A")

# ------------------------------------------------------------------------------
# (G) Removing and adding sidechains
Peptides.remove_sidechains!(pose, res_lib)
Peptides.add_sidechains!(pose, res_lib)

# Note: You can also specify a region/selection to remove/add the sidechains.
Peptides.remove_sidechains!(pose, res_lib, rid"4:5")
Peptides.add_sidechains!(pose, res_lib, rid"4:5")

# ------------------------------------------------------------------------------
# (H) Removing a residue in the middle of the segment
# Note: Since this manipulation breaks the segment, the downstream residues are
# not connected to the origin. Therefore, any changes upstream are effectively
# stopped at this point, not being passed downstream.
Peptides.pop_residue!(pose, pose.graph[1][4])

# ------------------------------------------------------------------------------
# (I) Setting a specific dihedral
# Note: `setdihedral!` function requires the requested angle to be in radians.
# By using ProtoSyn.Units module, several auxiliary tools are made available,
# such as the `°` operator, that facilitates the definition of radian angles.
# The Peptides.Dihedral module defines several function for the various angles
# in a peptide, receiving a Residue instance as input and outputting the
# atom responsible for setting that dihedral.
ProtoSyn.setdihedral!(pose.state, Peptides.Dihedral.phi(pose.graph[1][2]), 20°)