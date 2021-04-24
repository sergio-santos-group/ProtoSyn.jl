```@meta
CurrentModule = ProtoSyn
```

# [Pose](@id pose-methods)

As previously stated, a [Pose](@ref) is the main object of ProtoSyn and brings together both the [Graph](@ref graph-types) and the [State](@ref state-types). The following methods directly operate on [Pose](@ref) structs. These methods are subdivided by topics, for organization purposes:

+ [Creating, appending and inserting Fragments](@ref)
+ [Appending, inserting and removing Atom and Residue instances](@ref)
+ [Copying and merging](@ref)
+ [Aligning poses](@ref)

## Creating, appending and inserting Fragments

The next methods allow a user to create a new [Fragment](@ref) from a [Pose](@ref) instance and append or insert in in an existing structure. Note: creation of a new [Fragment](@ref) from a sequence of [`Residue`](@ref) instances is available at the [Builder](@ref) submodule.

```@docs
fragment(::Pose{Topology})
isfragment
append_fragment_as_new_segment!(::Pose{Topology}, ::Fragment)
append_fragment!(::Pose{Topology}, ::Residue, ::LGrammar, ::Pose{Segment}; ::Any)
insert_fragment!(::Pose{Topology}, ::Residue, ::LGrammar, ::Pose{Segment}; ::Any, ::Bool)
```

## Appending, inserting and removing Atom and Residue instances

In this section the methods responsible to manipulating the molecular structure of a [Pose](@ref) are explored.

```@docs
insert_atom_as_children!
pop_atom!
pop_residue!
```

## Copying and merging

The following methods allow a user to save or replicate an existing [Pose](@ref)
(for example, for simulations with multiple molecules).

```@docs
Base.copy(::Pose)
merge
merge!
symexp
symexp!
```

![ProtoSyn symexp](../../../assets/ProtoSyn-symexp.png)

**Figure 1 |** An example of the output of [`symexp`](@ref) with just one symmetry dimension. The second [Pose](@ref) is a copy of the original, translated in 1 dimension by the defined value.

## Aligning poses

The following methods allow a user to juxtapose two [`Pose`](@ref) instances and measure the RMSD distance between them.

```@docs
align!
center_of_mass
rmsd
```