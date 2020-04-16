export Opt
export Atom, Residue, Segment, Topology
# export AtomGraphNode, ResidueGraphNode

const Opt = Union{Nothing, T} where T

# should implemente the name and id fields
# abstract type AbstractIdentifiable end

# should implement the items[], container and size fields
abstract type AbstractContainer{T} end

abstract type AbstractAtom <: AbstractContainer{Nothing} end
abstract type AbstractResidue <: AbstractContainer{AbstractAtom} end
abstract type AbstractSegment <: AbstractContainer{AbstractResidue} end
abstract type AbstractTopology <: AbstractContainer{AbstractSegment} end


#-------------------------------------------------------------------------------
#region GraphNode
# mutable struct GraphNode{T}
#     item::T
#     parent::Opt{GraphNode{T}}
#     children::Vector{GraphNode{T}}
#     visited::Bool
#     ascendents::Opt{NTuple{4,Int}}
# end

# GraphNode{T}(item::T) where T = 
#     GraphNode(item, nothing, GraphNode{T}[], false, nothing)

#endregion GraphNode

function initgraph!(c::AbstractContainer)
    c.parent = nothing
    c.children = Atom[]
    c.visited = false
    c.ascendents = nothing
    c
end

#-------------------------------------------------------------------------------
#region Identifiables
mutable struct Atom <: AbstractAtom
    name::String                    # atom name
    id::Int                         # atom ID
    index::Int                      # global atom index (1-based)
    symbol::String                  # element symbol
    bonds::Vector{Atom}             # list of connected atoms
    container::Opt{AbstractResidue}   # parent residue
    # node::GraphNode{Atom}           #
    #------
    parent::Opt{Atom}
    children::Vector{Atom}
    visited::Bool
    ascendents::Opt{NTuple{4,Int}}
    #------
    Atom(name::String, id::Int, index::Int, symbol::String) = begin
        at = new(name, id, index, symbol, Atom[], nothing)
        # at.node = GraphNode{Atom}(at)
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
    # node::GraphNode{Residue}        #
    #------
    parent::Opt{Residue}
    children::Vector{Residue}
    visited::Bool
    ascendents::Opt{NTuple{4,Int}}
    #------
    Residue(name::String, id::Int) = begin
        r = new(name, id, Atom[], Dict{String,Atom}(), nothing, 0)
        # r.node = GraphNode{Residue}(r)
        initgraph!(r)
        r
    end
end

mutable struct Segment <: AbstractSegment
    name::String                    # segment name
    id::Int                         # segment ID
    items::Vector{Residue}       # list of residues (children)
    container::Opt{AbstractTopology} # parent topology
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
    Topology(name::String, id::Int) = begin
        new(name, id, Segment[], 0, Root())
    end
end

function Root()
    root = Residue("ROOT", -1)
    y = Atom("OY", -2, -2, "?")
    x = Atom("OX", -1, -1, "?")
    o = Atom("OO",  0,  0, "?")
    setparent!(y, x)
    setparent!(x, o)

    push!(root, y, x, o)
    root
end

#endregion Identifiable

# const AtomGraphNode = GraphNode{Atom}
# const ResidueGraphNode = GraphNode{Residue}

