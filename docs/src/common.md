# Common

```@meta
CurrentModule = Common
```

## Components

This section provides a description of the Common components.

### Energy

Contains common and simple energy representations. More specific energy structures can be used from other modules such as [`Amber`](@ref Forcefield)

```@docs
NullEnergy
Energy
```

### Residue

In ProtoSyn, a Residue object is a collection of atoms (normally an aminoacid) that are identified by a name and are part of a continuous tree of other Residues (have a `next` Residue).

```@docs
Residue
```

### Dihedral

A Dihedral is a collection of 4 atoms that define a dihedral in the simulated molecule. Many [Mutators](@ref Mutators) operate over this dihedrals, changing them in order to explore the conformational space of the system. A Dihedral is part of a [`Residue`](@ref) and has a defined DIHEDRAL.TYPE.

```@docs
Dihedral
rotate_dihedral!
rotate_dihedral_to!
```

### Metadata

ProtoSyn Metadata defines additional information of the system that is not necessarily requeired for the basic functions of the library, but allows for a better representation of the system. The [`AtomMetadata`](@ref) structure holds information related to each atom in the system, such as its element, connections, etc, while the [`SecondaryStructureMetadata`](@ref) holds informations regarding the system secondary structure organization.

```@docs
Metadata
AtomMetadata
SecondaryStructureMetadata
BlockMetadata
renumber_residues!
```

### State

The system state holds information about the current coordinates, energy and forces, aswell as any additional metadata. If iterated over, it returns atom by atom position and metadata.

```@docs
State
```

### Callback

The CallbackObject allows for independent calls to various functions with individual frequency of output. 

```@docs
CallbackObject
```

## Loaders

This section provides a description on how to load a new [`State`](@ref), [`Residue`](@ref) and [`Dihedral`](@ref) arrays, including the system [`Metadata`](@ref).

```@docs
load_from_gro
load_from_pdb
load_metadata_from_json
compile_residue_metadata!
compile_dihedral_metadata!
compile_ss
compile_blocks
```

## Conformation Generators

Conformation generators are responsible to change the system [`State`](@ref) in a limited and defined way.

```@docs
apply_ss!
apply_dihedrals_from_file!
stretch_conformation!
fix_proline!
```

## Macros

Auxiliary functions that help speed up the system's performance.

```@docs
@cbcall
@callback
@faggregator
```
