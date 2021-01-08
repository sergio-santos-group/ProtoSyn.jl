# In this example we will explore what types of selections ProtoSyn has
# available. Selections are auxiliary objects/macros that, when evaluated,
# return a subset of atoms/residues/segments that can then be used in multiple
# functions to restrict the application of certain manipulations to that subset.
# Note: For this example we will use the following built peptide as a canvas.
# ProtoSyn makes available the `print_selection` function, that helps in the
# visualization of the selected areas in a structure.

using ProtoSyn
using ProtoSyn.Peptides
using ProtoSyn.Builder

res_lib = Peptides.grammar()
pose = Peptides.build(res_lib, seq"AAAGGGKKKLLL")
ProtoSyn.write(pose, "selection.pdb")

# ------------------------------------------------------------------------------
# (A) By index and ID (and overall introduction to selections)
# Note: A selection is always an AbstractSelection object. This type is
# parametrized by the selection target (Atom, Residue, Segment, etc). This
# states the "level of detail" requested. SerialSelection is a type of selection
# that look for the `id` or `index` fields of the Atom, Residue, etc
selection = SerialSelection{Residue}(10, :index)
selection = SerialSelection{Atom}(10, :id)

# This selection can then be applied to a pose to get a Mask of the same type
mask = selection(pose)
ProtoSyn.print_selection(pose, mask, "selection.pdb")

# The mask can be used in the `ProtoSyn.gather` function to retrieve the list
# of item (Atoms, Residues, etc) that are set to `true` in the Mask.
atoms = ProtoSyn.gather(mask, pose.graph)

# These two steps can be condensed when applying the selection to a pose by
# setting the `gather` flag to `true`.
atoms = selection(pose, gather = true)

# Most selection types have short syntax macros. For example, for the two
# selections above:
selection = rid"10"
selection = aid"10"


# ------------------------------------------------------------------------------
# (B) Promotion

# One can also promote masks from any type to any other type using the `promote`
# function. This can be done downwards (Residue to Atom level, for example):
mask1 = rid"10"(pose) # Residue level
mask2 = ProtoSyn.promote(mask1, Atom, pose.graph) # Atom level

# Or upwards (Atom level to Residue level, for example). In this case, an
# aggregator function must be provided. The two most common aggregator functions
# are `any` (used by default, if any low-level object belonging to the
# high-level object is set to `true`) and `all` (all low-level objects must be
# `true` in the mask).
mask1 = aid"10"(pose) # Atom level
mask2 = ProtoSyn.promote(mask1, Residue, pose.graph, any) # Residue level

# Another useful promotion function returns two input masks promoted to the
# lowest ranking type. For example, if given an Atom and a Residue mask, the two
# output masks will be the Atom mask (untouched) and the Residue mask promoted
# to Atom level
mask1 = aid"10"(pose)
mask2 = rid"10"(pose)
new_mask1, new_mask2 = ProtoSyn.promote(mask1, mask2, pose.graph)

# The promotion step can also be pre-set on the selection itself, yielding a
# PromoteSelection type of selection. Same as before, an aggregator function
# should be specified in promoting upwards (`any` is used here, as default).
selection = ProtoSyn.promote(rid"10", Atom)


# ------------------------------------------------------------------------------
# (C) Binary Selections
# Selections can also be combined in `and` and `or` binary operations, using the
# `&` and `|` operations, respectively.
selection = BinarySelection(&, aid"10", aid"9")
selection = BinarySelection(|, aid"10", aid"9")

# Same as before, these operations have a short syntax version:
selection = aid"10" & aid"9"
selection = aid"10" | aid"9"

# When multiple combinations are done, ProtoSyn groups binary selection starting
# from the right, i.e: A | B | C => A | (B | C)
selection = aid"10" | aid"9" | aid"8"

# This behaviour can be overwritten be using parenthesis
selection = (aid"10" | aid"9") | aid"8"

# When combining selection of different levels (Atoms with Residues), the
# resulting selection is always of the lowest level of detail
selection = aid"10" & rid"1"


# ------------------------------------------------------------------------------
# (D) Unary Selection
# This is essentially the "negative" of a given selection.
selection1 = aid"10"
selection2 = !selection1
ProtoSyn.print_selection(pose, selection2(pose), "selection.pdb")


# ------------------------------------------------------------------------------
# (E) Field Selections
# These function search for specific fields in the Atom, Residue, etc objects.
# These are, essentially, the same as the previous SerialSelections, altough
# specialized for non-numeric entries, and as such also have short syntax
# versions.
selection = FieldSelection{Residue}("ALA", :name)
selection = rn"ALA"
ProtoSyn.print_selection(pose, selection(pose), "selection.pdb")

# These selection types can also use regular expressions (Regex) to search the
# molecular structure more efficiently, by setting the `is_regex` flag or adding
# the `r` flag to the short syntax.
selection = FieldSelection{Residue}("Y", :name, is_regex = true)
selection = rn"Y"r
ProtoSyn.print_selection(pose, selection(pose), "selection.pdb")

# As such, certain binary operations can be set inside the regular expression
# itself, such as:
selection = rn"GLY|ALA"r

# Other FieldSelection types include:
selection = sn"UNK" # By segment name
selection = an"UNK" # By atom name
selection = as"UNK" # By atom symbol


# ------------------------------------------------------------------------------
# (F) Distance Selections
# Distance selection select all Atom instances that are within X Angstroms of
# another selection (of Atoms, Residues, etc). Notice that the output mask is
# always of Atom type. It can later be promoted to any other type.
selection = DistanceSelection(10.0, rn"ALA")
ProtoSyn.print_selection(pose, selection(pose), "selection.pdb")

# The short syntax version is the following
selection = 10.0:rn"ALA"

# The selected part can be later masked out using the binary operations, for
# example. Notice the parenthesis on the distance selection. This is because, as
# previously stated, ProtoSyn employs right-first aggregation of selections.
# This means that, without the parenthesis, ProtoSyn would look for Atoms within
# 10 Angstrom of the 'rn"ALA" & !rn"ALA"' selection (which is obviously and
# empty set).
selection = (10.0:rn"ALA") & !rn"ALA"
ProtoSyn.print_selection(pose, selection(pose), "selection.pdb")

# ------------------------------------------------------------------------------
# (G) True Selections
# It can sometimes be helpful to have a selection/mask with all Atom, Residues,
# etc instances set to true.
selection = TrueSelection{Atom}()
selection = TrueSelection{Residue}()

# ------------------------------------------------------------------------------
# (H) Random Selections
# Random selections select a random Atom, Residue, etc from a container.
selection = RandomSelection{Residue}()
ProtoSyn.print_selection(pose, selection(pose), "selection.pdb")

# Alternatively, random selections can also take another selection from which to
# sample a random Atom, Residue, etc. In this example, we will select a random
# atom from an Alanine residue. 
selection = RandomSelection{Atom}(rn"ALA")
ProtoSyn.print_selection(pose, selection(pose), "selection.pdb")

# Another type of random selection is the RandomRangeSelection. This selects two
# random Atoms, Residues, etc from the given container and selects all instances
# between the two (based on the ID).
selection = RandomRangeSelection{Residue}()
ProtoSyn.print_selection(pose, selection(pose), "selection.pdb")