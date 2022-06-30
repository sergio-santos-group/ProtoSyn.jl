export Atom, Residue, Segment, Topology

"""
    AbstractContainer{T}

Supertype for a container of elements of type `T`.
"""
abstract type AbstractContainer{T} end

# abstract type AbstractAtom         <: AbstractContainer{Nothing}         end
# abstract type AbstractResidue      <: AbstractContainer{AbstractAtom}    end
# abstract type AbstractSegment      <: AbstractContainer{AbstractResidue} end
# abstract type AbstractTopology     <: AbstractContainer{AbstractSegment} end


function initgraph!(c::T) where {T <: AbstractContainer}
    c.children   = Vector{T}()
    c.ascendents = nothing
    c.parent     = nothing
    c.visited    = false
    c
end


"""
    Atom <: AbstractContainer{Nothing}

An `Atom` type.

    Atom(name::String, id::Int, index::Int, symbol::String)
    
Construct an `Atom`, with the given `name`, `id`, `index`, and `symbol`.
The created atom has no bonds and container; it is also a free directed-graph
node (no parent/children). *Note:* [`Atom`](@ref) is of super type
`AbstractAtom`, which is an `AbstractContainer`.

# Fields
* `name::String` - The name of the [`Atom`](@ref)
* `id::Int` - The ID of the [`Atom`](@ref)
* `index::Int` - The index of the [`Atom`](@ref) (default: same as `:id`)
* `symbol::String` - The chemical element symbol of this [`Atom`](@ref)
* `bonds::Vector{Atom}` - A list of [`Atom`](@ref) instances connected to this [`Atom`](@ref) by a bond (default: empty)
* `container::Opt{AbstractResidue}` - An optional container for this [`Atom`](@ref) (default: `nothing`)
* `visited::Bool` - Check whether this [`Atom`](@ref) has been visited (used by a some functions) (default: `false`)
* `parent::Opt{Atom}` - Optionally, the parent of this [`Atom`](@ref) in the directional graph (default: `nothing`)
* `children::Vector{Atom}` - Optionally, the list of children [`Atom`](@ref) instances downstream of this [`Atom`](@ref) (default: empty)
* `ascedents::Opt{NTuple{4, Int}}` - Optionally, the list of 4 ascendents, including this [`Atom`](@ref) (default: `nothing`)

# See also
[`Residue`](@ref) [ascendents]

# Examples
```jldoctest
julia> at = Atom("H1", 1, 1, "H")
Atom{/H1:1}
```
"""
mutable struct Atom <: AbstractContainer{Nothing}
    name::String                            # atom name
    id::Int                                 # atom ID
    index::Int                              # global atom index (1-based)
    symbol::String                          # element symbol
    bonds::Vector{Atom}                     # list of connected atoms
    container::Opt{AbstractContainer{Atom}} # parent residue
    
    visited::Bool
    parent::Opt{Atom}
    children::Vector{Atom}
    ascendents::Opt{NTuple{4,Int}}
    
    Atom(name::String, id::Int, index::Int, symbol::String) = begin
        at = new(name, id, index, symbol, Atom[], nothing)
        initgraph!(at)
        at
    end
end


"""
    Residue <: AbstractContainer{Atom}
    Residue(name::String, id::Int)
    
Construct a [`Residue`](@ref), with the given `name` and `id`. The created
residue has no container; it is also a free directed-graph node
(no parent/children). A [`Residue`](@ref) is a group of [`Atom`](@ref)
instances. *Note:* [`Residue`](@ref) is of super type `AbstractResidue`, which
is an `AbstractContainer`.

# Fields
* `name::ResidueName` - The name of the [`Residue`](@ref)
* `id::Int` - The ID of the [`Residue`](@ref)
* `index::Int` - The index of the [`Residue`](@ref) (default: same as `:id`)
* `items::Vector{Atom}` - A list of [`Atom`](@ref) instances in this [`Residue`](@ref) (default: empty)
* `itemsbyname::Dict{String, Atom}` - A `Dict` object containing a link to each of the [`Atom`](@ref) instances in this [`Residue`](@ref) bsed on the `:name` of the [`Atom`](@ref) (default: empty)
* `container::Opt{AbstractSegment}` - An optional container for this [`Residue`](@ref) (default: `nothing`)
* `size::Int` - The number of [`Atom`](@ref) instances in this [`Residue`](@ref) (default: 0)
* `visited::Bool` - Check whether this [`Residue`](@ref) has been visited (used by a some functions) (default: `false`)
* `parent::Opt{Residue}` - Optionally, the parent of this [`Residue`](@ref) in the directional graph (default: `nothing`)
* `children::Vector{Residue}` - Optionally, the list of children [`Residue`](@ref) instances downstream of this [`Residue`](@ref) (default: empty)
* `ascedents::Opt{NTuple{4, Int}}` - Optionally, the list of 4 [`ascendents`](@ref), including this [`Residue`](@ref) (default: `nothing`)

# See also
[`Segment`](@ref) [`Atom`](@ref) [`ascendents`](@ref) [`ResidueName`](@ref)

# Examples
```jldoctest
julia> res = Residue("UNK", 1)
Residue{/UNK:1}
```
"""
mutable struct Residue <: AbstractContainer{Atom}
    name::ResidueName                          # residue name
    id::Int                                    # residue ID
    index::Int                                 # residue index
    items::Vector{Atom}                        # list of atoms
    itemsbyname::Dict{String, Atom}            # list of atoms indexed by name
    container::Opt{AbstractContainer{Residue}} # container segment
    size::Int
    
    visited::Bool
    parent::Opt{Residue}
    children::Vector{Residue}
    ascendents::Opt{NTuple{4,Int}}
    
    Residue(name::String, id::Int) = begin
        r = new(ResidueName(name), id, id, Atom[], Dict{String,Atom}(), nothing, 0)
        initgraph!(r)
        r
    end
end


"""
    Segment <: AbstractContainer{Residue}
    Segment(name::String, id::Int)
    
Construct a [`Segment`](@ref), with the given `name` and `id`. The created
segment has no container. A [`Segment`](@ref) is a group of [`Residue`](@ref)
instances. *Note:* [`Segment`](@ref) is of super type `AbstractSegment`, which
is an `AbstractContainer`.

# Fields
* `name::String` - The name of the [`Segment`](@ref)
* `id::Int` - The ID of the [`Segment`](@ref)
* `index::Int` - The index of the [`Segment`](@ref) (default: same as `:id`)
* `code::Char` - A 1-letter code for the chain (used in PDB formatting) (default: '?')
* `items::Vector{Residue}` - A list of [`Residue`](@ref) instances in this [`Segment`](@ref) (default: empty)
* `container::Opt{AbstractTopology}` - An optional container for this [`Segment`](@ref) (default: `nothing`)
* `size::Int` - The number of [`Residue`](@ref) instances in this [`Segment`](@ref) (default: 0)

# See also
[`Topology`](@ref) [`Residue`](@ref)

# Examples
```jldoctest
julia> res = Segment("UNK", 1)
Segment{/UNK:1}
```
"""
mutable struct Segment <: AbstractContainer{Residue}
    name::String                               # segment name
    id::Int                                    # segment ID
    index::Int                                 # segment index
    code::Char
    items::Vector{Residue}                     # list of residues
    container::Opt{AbstractContainer{Segment}} # container topology
    size::Int
    
    Segment(name::String, id::Int) = begin
        new(name, id, id, '?', Residue[], nothing, 0)
    end
end


"""

    Topology <: AbstractContainer{Segment}
    Topology(name::String, id::Int)
    
Construct a [`Topology`](@ref), with the given `name` and `id`. A
[`Topology`](@ref) is the top level of hierarchy in a Graph, and holds a `Root`
[`Residue`](@ref), therefore initializing the internal coordinates system in the
corresponding [State] of a [Pose](@ref pose-types).  *Note:* [`Topology`](@ref) is of super
type `AbstractTopology`, which is an `AbstractContainer`.

# Fields
* `name::String` - The name of the [`Topology`](@ref)
* `id::Int` - The ID of the [`Topology`](@ref)
* `items::Vector{Segment}` - All the [`Segment`](@ref) instances in this [`Topology`](@ref) (default: empty)
* `size::Int` - The number of [`Segment`](@ref) instances in this [`Topology`](@ref) (default: 0)
* `root::Residue` - The root [`Residue`](@ref) of this [`Topology`](@ref)
* `container::Nothing` - Used to identify the root [`Residue`](@ref)  (default: `nothing`)

!!! ukw "Note:"
    When in a [Pose](@ref pose-types), both the [`Topology`](@ref) and the corresponding
    [State](@ref state-types) need to share the same `:id`. This is used to identify
    situations where a change in one of this structures was not reflected on the
    other.

# See also
[`root`](@ref) [`Segment`](@ref)

# Examples
```jldoctest
julia> top = Topology("UNK", 1)
Topology{/UNK:1}
```
"""
mutable struct Topology <: AbstractContainer{Segment}
    name::String
    id::Int
    items::Vector{Segment}
    size::Int
    root::Residue
    
    container::Nothing
    
    Topology(name::String, id::Int) = begin
        new(name, id, Segment[], 0, Root(), nothing)
    end
end


#region order ----------------------------------------------------------------

# AbstractContainer order
# Allows for "Atom < Residue" and "min(Atom, Residue)", as examples
Base.isless(::Type{<: T}, ::Type{<: AbstractContainer{T}}) where {T <: AbstractContainer} = true
Base.isless(::Type{<: T}, ::Type{<: AbstractContainer{<: AbstractContainer{T}}}) where {T <: AbstractContainer} = true
Base.isless(::Type{<: AbstractContainer{T}}, ::Type{<: T}) where {T <: AbstractContainer} = false
Base.isless(::Type{<: AbstractContainer{<: AbstractContainer{T}}}, ::Type{<: T}) where {T <: AbstractContainer} = false
Base.isless(::Type{T}, ::Type{T}) where {T <: AbstractContainer} = false
Base.isless(::Type{Atom}, ::Type{Topology}) = true
Base.isless(::Type{Topology}, ::Type{Atom}) = false

# endregion order


"""
    Root()::Residue
    
Return a new Root residue. A Root residue in a pseudo-residue that serves the
purpose of establishing an anchor for initial internal coordinates definition.

# Examples
```jldoctest
julia> root = ProtoSyn.Root()
Residue{/ROOT:17124}
```
"""
function Root()::Residue
    root = Residue("ROOT", genid())
    y = Atom("OY", -2, -2, "?")
    x = Atom("OX", -1, -1, "?")
    o = Atom("OO",  0,  0, "?")
    push!(root, y, x, o)
    setparent!(x, y)
    setparent!(o, x)
    root
end


export Segment!, Residue!, Atom!

"""
    Segment!(topology::Topology, name::String, id::Int)::Residue
    
Create a new `Segment` (with the given `name` and `id`) and add it to the given
Topology `topology`. Returns the created `Segment` instance.

# Examples
```jldoctest
julia> seg = Segment!(pose.graph, "UNK", 1)
Segment{/UNK:1/UNK:1}
```
"""
function Segment!(topology::Topology, name::String, id::Int)::Segment
    s = Segment(name, id)
    push!(topology, s)
    s
end


"""
    Residue!(segment::Segment, name::String, id::Int)::Residue
    
Create a new `Residue` (with the given `name` and `id`) and add it to the given
Segment `segment`. Returns the created `Residue` instance.

# Examples
```jldoctest
julia> res = Residue!(pose.graph[1], "ALA", 1)
Residue{/UNK:1/UNK:1/ALA:1}
```
"""
function Residue!(segment::Segment, name::String, id::Int)::Residue
    r = Residue(name, id)
    push!(segment, r)
    r
end

Residue!(segment::Segment, name::ResidueName, id::Int) = begin
    Residue!(segment, name.content, id)
end


"""
    Atom!(r::Residue, name::String, id::Int, index::Int, symbol::String)
    
Create a new `Atom` (with the given `name`, `id`, `index` and `symbol`) and add
it to the given Residue `residue`. Returns the created `Atom` instance.

# Examples
```jldoctest
julia> atom = Atom!(pose.graph[1][1], "H1", 1, 1, "H")
Atom{/UNK:1/UNK:1/SER:1/H1:1}
```
"""
function Atom!(r::Residue, name::String, id::Int, index::Int, symbol::String)::Atom
    a = Atom(name, id, index, symbol)
    push!(r, a)
    a
end