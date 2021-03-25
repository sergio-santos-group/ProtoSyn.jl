export Opt
export Atom, Residue, Segment, Topology

const Opt = Union{Nothing, T} where T

"""
    AbstractDigraph

Supertype for a directed graph node. Types implementing the AbstractDigraph
interface should have the following fields:
 - parent::Union{T,Nothing}
 - children::Vector{T}
 - visited::Bool
 - ascendents::Union{Nothing,NTuple{N,Int}}
"""
abstract type AbstractDigraph end


"""
    AbstractContainer{T} <: AbstractDigraph

Supertype for a container of elements of type `T`, implementing
the `AbstractDigraph` interface.
"""
abstract type AbstractContainer{T} <: AbstractDigraph  end

abstract type AbstractAtom         <: AbstractContainer{Nothing}         end
abstract type AbstractResidue      <: AbstractContainer{AbstractAtom}    end
abstract type AbstractSegment      <: AbstractContainer{AbstractResidue} end
abstract type AbstractTopology     <: AbstractContainer{AbstractSegment} end


function initgraph!(c::T) where {T<:AbstractDigraph}
    c.children   = Vector{T}()
    c.ascendents = nothing
    c.parent     = nothing
    c.visited    = false
    c
end


"""
    Atom <: AbstractAtom

An `Atom` type.

    Atom(name::String, id::Int, index::Int, symbol::String)
    
Construct an `Atom`, with the given `name`, `id`, `index`, and `symbol`.
The created atom has no bonds and container; it is also a free directed-graph
node.

# Examples
```jldoctest
julia> at = Atom("H1", 1, 1, "H")
Atom{/H1:1}

julia> hascontainer(at), isempty(at.bonds)
(false, true)

julia> hasparent(at), haschildren(at)
(false, false)
```
"""
mutable struct Atom <: AbstractAtom
    name::String                    # atom name
    id::Int                         # atom ID
    index::Int                      # global atom index (1-based)
    symbol::String                  # element symbol
    bonds::Vector{Atom}             # list of connected atoms
    container::Opt{AbstractResidue} # parent residue
    
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

segment(at::Atom) = hascontainer(at) ? at.container.container : nothing


"""
    Residue <: AbstractResidue

A `Residue` type.

    Residue(name::String, id::Int)
    
Construct a `Residue`, with the given `name` and `id`.
The created residue has no container; it is also a free directed-graph
node.

# Examples
```jldoctest
julia> res = Residue("UNK", 1)
Residue{/UNK:1}
```
"""
mutable struct Residue <: AbstractResidue
    name::ResidueName                    # residue name
    id::Int                         # residue ID
    index::Int                      # residue index
    items::Vector{Atom}             # list of atoms (children)
    itemsbyname::Dict{String, Atom} # child atoms indexed by name
    container::Opt{AbstractSegment} # parent segment
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
    Segment <: AbstractSegment

A `Segment` type.

    Segment(name::String, id::Int)
    
Construct a `Segment`, with the given `name` and `id`.
The created segment has no container.

# Examples
```jldoctest
julia> res = Segment("UNK", 1)
Segment{/UNK:1}
```
"""
mutable struct Segment <: AbstractSegment
    name::String                        # segment name
    id::Int                             # segment ID
    index::Int                          # segment index
    code::Char
    items::Vector{Residue}              # list of residues (children)
    container::Opt{AbstractTopology}    # parent topology
    size::Int
    
    Segment(name::String, id::Int) = begin
        new(name, id, id, '?', Residue[], nothing, 0)
    end
end


"""
    Topology <: AbstractTopology

A `Topology` type.

    Topology(name::String, id::Int)
    
Construct a `Topology`, with the given `name` and `id`.

# Examples
```jldoctest
julia> res = Segment("UNK", 1)
Segment{/UNK:1}
```
"""
mutable struct Topology <: AbstractTopology
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
Residue{/ROOT:53033}
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
julia> seg = Segment!(topology, "UNK", 1)
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
julia> res = Residue!(segment, "ALA", 1)
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
julia> atom = Atom!(residue, "H1", 1, 1, "H")
Atom{/UNK:1/UNK:1/ALA:1/H1:1}
```
"""
function Atom!(r::Residue, name::String, id::Int, index::Int, symbol::String)::Atom
    a = Atom(name, id, index, symbol)
    push!(r, a)
    a
end