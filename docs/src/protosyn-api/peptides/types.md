```@meta
CurrentModule = ProtoSyn
```


# [Types] (@id peptides-types)

The [Peptides](@ref) module introduces two main types: the [Dihedral](@ref) and the [Secondary Structure](@ref). In the next section, a more in depth exploration of the usage and utility of each of these types is provided, divided in the following subsections for organization purposes:

+ [Dihedral](@ref)
+ [Secondary Structure](@ref)

# Dihedral

The [Dihedral](@ref) struct allows a user to select a given dihedral in a [Pose](@ref) by it's name in the context of a protein structure (i.e.: select the "phi" dihedral or the "chi-2" dihedral). This is programmed as a Julia module for sugar syntax.

```@docs
ProtoSyn.Peptides.Dihedral
```

```@meta
CurrentModule = ProtoSyn.Peptides
```

![ProtoSyn Chi Dihedral Atoms](../../assets/ProtoSyn-chi-dihedral-atoms.png)

**Figure 1 |** In ProtoSyn, peptide atoms are named according to the IUPAC and Protein Data Bank nomeculature (see the [Atomic Coordinate and Bibliographic Entry Format Description](https://cdn.rcsb.org/wwpdb/docs/documentation/file-format/PDB_format_1992.pdf) - page 26). As such, in the table above, the names of the representative atoms of each chi angle are shown. These are the atoms returned by the call to the `Dihedral` auxiliary struct, and can be used in other functions, such as [`getdihedral`](@ref ProtoSyn.getdihedral) or [`setdihedral!`](@ref ProtoSyn.setdihedral!).

# Secondary Structure

The following types allow the maintenance of [Secondary Structure](@ref) information (i.e.: the `psi`, `psi` and `omega` backbone angles that result in any given [Secondary Structure](@ref), such as an alpha helix or a beta sheet, for example).

```@docs
SecondaryStructureTemplate
SecondaryStructure
```

![ProtoSyn Backbone SS](../../assets/ProtoSyn-backbone-ss.png)

**Figure 2 |** The default [`SecondaryStructureTemplate`](@ref) instances in ProtoSyn, defined in the [`SecondaryStructure`](@ref) dictionary. A [`SecondaryStructureTemplate`](@ref) can be applied to a [`Pose`](@ref), for example, using the [`setss!`](@ref) method.