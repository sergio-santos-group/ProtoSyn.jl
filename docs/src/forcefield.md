# Forcefield

Currently, ProtoSyn only supports the Amber forcefield.

```@meta
CurrentModule = Forcefield
```

## Components

This section provides a description on all the individual components that comprise the Forcefield.


```@docs
Amber.HarmonicBond
Amber.HarmonicAngle
Amber.DihedralCos
Amber.Atom
Amber.Topology
```

## Evaluators

This section provides a description on how to use [Forcefield.Components](#Components-1) to calculate their energy.


```@docs
Amber.evaluate!
```

## Loaders

This section provides a description on how to load [Forcefield.Components](#Components-1) from external sources.

```@docs
Amber.load_from_json
```

