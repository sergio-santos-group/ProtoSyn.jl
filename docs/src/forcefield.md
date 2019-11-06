# Forcefield

```@meta
CurrentModule = Forcefield
```

## Amber

### Components

This section provides a description on all the individual components that comprise the Amber Forcefield.


```@docs
Amber.HarmonicBond
Amber.HarmonicAngle
Amber.DihedralCos
Amber.Atom
Amber.Topology
```

### Evaluators

This section provides a description on how to use [Forcefield.Components](#Components-1) to calculate their energy.


```@docs
Amber.evaluate!
```

### Loaders

This section provides a description on how to load [Forcefield.Components](#Components-1) from external sources.

```@docs
Amber.load_from_json
```

## Restraints

### Components

This section provides a description on all the individual components that comprise the Restraints.


```@docs
Restraints.DistanceFBR
Restraints.DihedralFBR
```

### Evaluators

This section provides a description on how to use [Forcefield.Components](#Components-2) to calculate their energy.


```@docs
Restraints.evaluate!
```

### Loaders

This section provides a description on how to load [Forcefield.Components](#Components-2) from external sources.

```@docs
Restraints.load_distance_restraints_from_file
Restraints.lock_block_bb
```

## Coarse Grain

### Components

This section provides a description on all the individual components that comprise the Coarse Grain energy functions.


```@docs
CoarseGrain.SolvPair
CoarseGrain.HbGroup
```

### Evaluators

This section provides a description on how to use [Forcefield.Components](#Components-3) to calculate their energy.


```@docs
CoarseGrain.evaluate!
```

### Loaders

This section provides a description on how to load [Forcefield.Components](#Components-3) from external sources.

```@docs
CoarseGrain.load_solv_pairs_from_file
CoarseGrain.compile_solv_pairs
CoarseGrain.compile_hb_groups
```
