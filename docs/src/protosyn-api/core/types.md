# Types

In this section we will explore the main `struct` instances that compose the core of ProtoSyn's engine.

## Pose

A pose holds all necessary information regarding a molecular system. This information is divided in the [Graph](@ref) and [State](@ref) .
```@docs
Pose
```

## Fragment

ProtoSyn.jl package uses a pose of a single [`Segment`](@ref) (therefore without a root) to carry information during certain processes (such as [Building a new peptide from a sequence]).
```@docs
Fragment
```

## Graph

The birectional graph holds the information regarding the hierarchical organization of a molecular system, where the top level of organization is a [`Topology`](@ref) instance, which can contain one or more instances of [`Segment`](@ref). This structure, in turn, is comprised of one or more [`Residue`](@ref) instances, which are themselves a set of [`Atom`](@ref) instances, this being the lowest level of data organization in ProtoSyn.

![graph](../../assets/ProtoSyn-graph.png)

```@docs
Topology
Segment
Residue
```

!!! ukw "Note:"
    Instead of defining [`Residue`](@ref).name as a String, ProtoSyn uses a specific type ([`ResidueName`](@ref ProtoSyn.ResidueName)) to accomodate specific conventions in aminoacid naming.
```@docs
ProtoSyn.ResidueName
```

```@docs
Atom
```


## State

A [State](@ref) structure is responsible for holding the current _state_ of the internal and cartesian coordinates, as well as the forces felt on each atom and the resulting energy of the system (as calculated by an [EnergyFunction]). This structure is sub-divided in a list of [`AtomState`](@ref) instances and a [`StateMatrix`](@ref). Both of these structures are complementary, meaning that any change in one in synched in the other. While [`AtomState`](@ref) instances hold both the internal and cartesian coordinates relative to a single [`Atom`](@ref), [`StateMatrix`](@ref) is a 2D matrix of all cartesian coordinates of the atoms in a given [Pose](@ref).

```@docs
State
AtomState
StateMatrix
```