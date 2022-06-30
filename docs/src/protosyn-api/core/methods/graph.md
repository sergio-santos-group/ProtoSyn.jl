```@meta
CurrentModule = ProtoSyn
```

# [Graph](@id core-graph-methods)

This section lists functions that work on the system [Graph](@ref state-types), such as functions that deal with parent/child relations, indexes and bonds, among
others. These are subdivided by topics, for organization purposes:

+ [Root vs Origin](@ref)
+ [Parenthood relationships](@ref)
+ [Container manipulation](@ref)
+ [Indexation](@ref core-graph-methods-indexation)
+ [Counters and Iterators](@ref)
+ [Bonds](@ref)
+ [Travelling the Graph](@ref)

## Root vs Origin

The __root__ of a [Graph](@ref state-types) is a set of 3 pseudoatoms belonging to a [`Topology`](@ref), acting as the anchor for internal coordinate calculations, while an __origin__ is any [`Atom`](@ref) instance with no parent in any given container (for example, due to a cut). While a [`Pose`](@ref) can only have 1 __root__, it can have multiple __origins__ (which usually need to be reconnected for most of the simulations of ProtoSyn).

```@docs
origin
root
Root
```

## Parenthood relationships

[`Atom`](@ref) and [`Residue`](@ref) instances have parenthood relationships, in a directional [Graph](@ref state-types). Several methods allow the manipulation of such relationships, and are of additional importance when using internal coordinates to facilitate and speed up some calculations.

```@docs
hasparent
isparent
setparent!
popparent!
ascendents
Base.detach
is_contiguous
infer_parenthood!
```

![ProtoSyn graph](../../../assets/ProtoSyn-graph2.png)

**Figure 1 |** A diagram of a directional [Graph](@ref state-types) in ProtoSyn.

## Container manipulation

[`Atom`](@ref) instances have a `:container` field, setting the container `AbstractContainer` (usually a [`Residue`](@ref) instance). The same logic applies to [`Residue`](@ref) instances, inside [`Segment`](@ref) structs, and [`Segment`](@ref) instances inside [`Topology`](@ref) structs (according to the established hierarchy of `AbstractContainer` instances). The following methods allow for the correct manipulation of this relations, allowing to add and remove `AbstractContainer` instances from the respective `:container` field, as well as creating copies of them.

```@docs
hascontainer
Base.push!(::AbstractContainer{T}, ::T) where {T<:AbstractContainer}
Base.insert!(::AbstractContainer{T}, ::Integer, ::T) where {T <: AbstractContainer}
Base.delete!(::AbstractContainer{T}, ::T) where {T<:AbstractContainer}
Base.copy(::Atom)
rename!
Atom!
Residue!
Segment!
```

## [Indexation](@id core-graph-methods-indexation)

An important initial detail when describing the [Graph](@ref core-graph-methods) methods is describing the family of `getindex` methods overloaded by ProtoSyn when dealing with [Graph](@ref graph-types) structures. There is, in essence, 4 ways to access a specific instance in the [Graph](@ref graph-types):

* Using the regular syntax;

```@setup methods
using ProtoSyn
using ProtoSyn.Peptides
pose = ProtoSyn.build(Peptides.grammar, seq"GME")
```

```@repl methods
pose.graph.items[1].items[1].items[1]
pose.graph.items[1].items[1]
```

* Using the short syntax with index;

```@repl methods
pose.graph[1][1][1]
pose.graph[1][1]
```

* Using the condensed syntax with index;

```@repl methods
pose.graph[1, 1, 1]
pose.graph[1, 1]
```

* Using the short or condensed syntax with atom name (for [`Atom`](@ref) instances only).

```@repl methods
pose.graph[1][1]["N"]
pose.graph[1, 1, "N"]
```

Note that queries by [`Atom`](@ref)`.name` are case sensitive. Besides using it to get [Graph](@ref graph-types) structures, [`Atom`](@ref)`.name` fields can also be used to query for existence:

```@repl methods
"N" in pose.graph[1][1]
```

The following methods deal with the correct indexation of the [Graph](@ref state-types). Note that, altough not necessary, some simulation functions assume that both the [Graph](@ref state-types) and [State](@ref state-types) indexation of a [Pose](@ref pose-types) are synched and are equal.

```@docs
genid
reindex(::Topology; ::Bool)
ids
```

## Counters and Iterators

ProtoSyn includes custom [Counters and Iterators](@ref) to analyze and loop over [Graph](@ref graph-types) structures. In one hand, by using Counters the user can count the number of sub `AbstractContainer` instances in a [Graph](@ref graph-types) component (for example, the number of [`Atom`](@ref) instances in a [`Segment`](@ref) or a [`Residue`](@ref)). This can be achieved by using the `count_atoms`, `count_residues` and `count_segments` methods.

```@repl methods
ProtoSyn.count_atoms(pose.graph)
ProtoSyn.count_atoms(pose.graph[1][1])
ProtoSyn.count_residues(pose.graph)
```

On the other hand, Iterators allow the user to iterate over all sub `AbstractContainer` instances in a [Graph](@ref graph-types) component, by using the `eachatom`, `eachresidue` and `eachsegment` methods.

```@repl methods
eachatom(pose.graph)
for residue in eachresidue(pose.graph[1])
    println(residue);
end
```

For [`Residue`](@ref) instances in specific, a more direct way to list all instances in a given [`Pose`](@ref) or `AbstractContainer` is to use the [`sequence`](@ref) method.

```@docs
sequence
```

## Bonds

The following methods deal with the bonding/unbonding of atoms (and respective
[Parenthood relationships](@ref) when using the [`join`](@ref) function).

```@docs
bond
unbond!
join
infer_bonds!
```

## Travelling the Graph

As further explored in the [Graph](@ref graph-types) section, the directed nature of the [`Pose`](@ref)'s [Graph](@ref graph-types) allows for easy travelling of the system. The following methods facilitate that process.

```@docs
get_graph_size
sort_children
travel_graph
travel_bonds
identify_atom_by_bonding_pattern
```

As explained above, certain methods in ProtoSyn travel the directed graph. There are, however, multiple ways to follow the same graph (as longs as its ramified), depending on the criteria used on bifurcations. In the Core module, ProtoSyn makes available for the BFS (breath-first search) and DFs (depth-first search) algorithms, a type of `SearchAlgorithm`.

```@docs
BFS
DFS
```