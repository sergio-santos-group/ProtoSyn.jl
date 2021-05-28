```@meta
CurrentModule = ProtoSyn.Peptides
```

# [Pose](@id peptides-pose-methods)

Included in the [Peptides](@ref) are a set of methods who primarily act on a [`Pose`](@ref) instance, expanding on the list of availabe Core functions (See [Pose](@ref pose-methods) methods). These are subdivided by topics, for organization purposes:

+ [Appending and inserting fragments](@ref)
+ [Removing Atom and Residue instances](@ref)
+ [Mutating an aminoacid](@ref)
+ [Removing and adding sidechains](@ref)
+ [Removing and adding N- and C- terminals](@ref)

# Appending and inserting fragments

Expanding on the homologous Core functions, the next methods allow the user to append and insert [`Fragment`](@ref) instances in [`Pose`](@ref) instances, adding the possibility to set the [Secondary Structure](@ref) of the appendage (using the [`setss!`](@ref) method).

```@docs
append_fragment!
insert_fragment!
```

WIP

# Removing Atom and Residue instances

WIP

# Mutating an aminoacid

WIP

# Removing and adding sidechains

WIP

# Removing and adding N- and C- terminals

WIP