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

## Crankshaft

```@meta
CurrentModule = Mutators.Crankshaft
```

A Crankshaft movement is a rotation of all atoms between two randomly chosen alpha carbons according to a set of [`CrankshaftMutator`](@ref) parameters.

```@docs
CrankshaftMutator
run!
rotate_crankshaft!
```