```@meta
CurrentModule = ProtoSyn.Peptides
```

# [Builder](@id peptides-builder)

> The [Builder](@ref peptides-builder) is a submodule of `ProtoSyn.Peptides` module. As such, the following section introduces both new types and methods that work together, in a generally independent way from the rest of the module, and require an unique exploratory section on their own.

As an expansion of the Core module [Builder](@ref core-builder), this submodule introduces the peptidic [`LGrammar`](@ref) type, as well as the necessary methods to append and insert [`Fragment`](@ref) instances from a derivation, organized in the following topics:

+ [Loading the default Peptides Stochastic L-Grammar](@ref peptides-builder-1)
+ [Building a molecular structure](@ref peptides-builder-2)
+ [Manipulating a molecular structure by adding new residues from templates](@ref peptides-builder-3)

# [Loading the default Peptides Stochastic L-Grammar](@id peptides-builder-1)

```@docs
grammar(::Type{T}; ::Bool) where {T <: AbstractFloat}
```

# [Building a molecular structure](@id peptides-builder-2)

```@docs
build(::LGrammar{T}, ::Any, ::SecondaryStructureTemplate) where {T <: AbstractFloat}
```

# [Manipulating a molecular structure by adding new residues from templates](@id peptides-builder-3)

```@docs
append_fragment!(::Pose{Topology}, ::Residue, ::LGrammar, ::Any; ::Opt{SecondaryStructureTemplate}, ::Any)
insert_fragment!(::Pose{Topology}, ::Residue, ::LGrammar, ::Any; ::Opt{SecondaryStructureTemplate}, ::Any)
```