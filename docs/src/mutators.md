# Mutators

Mutators are modules that introduce a single conformational change on the system, according to a set of
runtime parameters. These can be combined into multiple steps, defining a [Driver](@ref Drivers).

## Dihedral

```@meta
CurrentModule = Mutators.Dihedral
```

This section provides a description on the Dihedral mutator, responsible for performing a single dihedral
movement according to a [`DihedralMutator`](@ref) set of parameters.

```@docs
DihedralMutator
run!
```