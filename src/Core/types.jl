export Opt
export Atom, Residue, Segment, Topology


const Opt = Union{Nothing, T} where T

# Subtypes of AbstractContainer shoud implement the following fields:
#  - name
#  - id
#  - items[]
#  - container
#  - size
abstract type AbstractContainer{T} <: AbstractDigraph  end
abstract type AbstractAtom         <: AbstractContainer{Nothing}         end
abstract type AbstractResidue      <: AbstractContainer{AbstractAtom}    end
abstract type AbstractSegment      <: AbstractContainer{AbstractResidue} end
abstract type AbstractTopology     <: AbstractContainer{AbstractSegment} end


function initgraph!(c::T) where {T<:AbstractDigraph}
    c.children = Vector{T}()
    c.ascendents = nothing
    c.parent = nothing
    c.visited = false
    c
end

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


mutable struct Residue <: AbstractResidue
    name::String                    # residue name
    id::Int                         # residue ID
    items::Vector{Atom}             # list of atoms (children)
    itemsbyname::Dict{String, Atom} # child atoms indexed by name
    container::Opt{AbstractSegment} # parent segment
    size::Int
    
    visited::Bool
    parent::Opt{Residue}
    children::Vector{Residue}
    ascendents::Opt{NTuple{4,Int}}
    
    Residue(name::String, id::Int) = begin
        r = new(name, id, Atom[], Dict{String,Atom}(), nothing, 0)
        initgraph!(r)
        r
    end
end


mutable struct Segment <: AbstractSegment
    name::String                        # segment name
    id::Int                             # segment ID
    items::Vector{Residue}              # list of residues (children)
    container::Opt{AbstractTopology}    # parent topology
    size::Int
    
    Segment(name::String, id::Int) = begin
        new(name, id, Residue[], nothing, 0)
    end
end


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

export ReactionToolbelt
struct ReactionToolbelt{F<:Function, G<:Function, H<:Function}
    join::F
    split::G
    root::H
end

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

