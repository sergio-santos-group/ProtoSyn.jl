```@meta
CurrentModule = ProtoSyn.Peptides
```

# [State](@id peptides-state-methods)

This section lists functions that work on the system [State](@ref state-types),
as expansions to the Core methods (See [State](@ref core-state-methods)). These methods are subdivided by topics, for organization purposes:

+ [Applying a Secondary Structure](@ref)

# Applying a Secondary Structure

Making use of the [Secondary Structure](@ref) types (in specific, the [`SecondaryStructureTemplate`](@ref) instance), the [`setss!`](@ref) method allows a user to quickly apply a given [Secondary Structure](@ref) to a [`Pose`](@ref) or a subset of selected [`Residue`](@ref) instances.

```@docs
setss!
```