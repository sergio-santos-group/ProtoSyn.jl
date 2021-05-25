# Peptides

Besides the Core module, ProtoSyn comes with several other extra modules that build on top of and expand the base methods, types and submodules. These ProtoSyn modules are specific for a given topic. In the case of the [Peptides](@ref) module, all new or expanded methods and types are reffering to the peptide and protein structures. With this extra constraint, most methods can be especialized: the [`ProtoSyn.bond`](@ref ProtoSyn.bond), for example, is specialized in [`ProtoSyn.Peptides.bond`](@ref ProtoSyn.Peptides.bond), automatically connecting the `C` atom to the next aminoacid's `N` atom - the peptidic bond.

In sum, the [Peptides](@ref) module makes available several methods, types and submodules specific for proteins and peptide structures.

In order to use this module, including the following call is often useful:

```@setup peptides
using ProtoSyn
```

```@example peptides
using ProtoSyn.Peptides
```