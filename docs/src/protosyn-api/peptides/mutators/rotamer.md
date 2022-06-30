```@meta
CurrentModule = ProtoSyn.Peptides.Mutators
```

# Rotamer Mutator

One of the `AbstractMutator` types made available by the ProtoSyn.Peptides module is the [`RotamerMutator`](@ref), allowing a user to sample the sidechain conformation of one or more [`Residue`](@ref) instances. This has, for example, applications in sidechain packaging, ligand adsorption stabilization or active site generation algorithms, among others. For more information, see the [Rotamers](@ref) submodule section.

```@docs
RotamerMutator
```

![ProtoSyn Rotamer Mutator Card](../../../assets/ProtoSyn-rotamer-mutator-card.png)

**Figure 1 |** Schematic description of the [`RotamerMutator`](@ref). 

![ProtoSyn Rotamer Mutator](../../../assets/ProtoSyn-rotamer-mutator.gif)

**Figure 2 |** Example of the application of [`RotamerMutator`](@ref) to simulateously change the conformation of all sidechain [`Atom`](@ref) instances (except in Proline [`Residue`](@ref) instances). 