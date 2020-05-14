export hasparent, haschildren, isparent
export setparent!, popparent!
export popchild!

# Types implementing the AbstractDigraph interface
# should have the following fields
#  - parent::Union{T,Nothing}
#  - children::Vector{T}
#  - visited::Bool
#  - ascendents::Union{Nothing,NTuple{N,Int}}
"""
    AbstractDigraph

Supertype for a directed graph node.
"""
abstract type AbstractDigraph end

"""
    hasparent(c::AbstractDigraph) -> Bool

Test whether the given AbstractDigraph `c` has a parent.
"""
@inline hasparent(c::AbstractDigraph) = c.parent !== nothing


"""
    haschildren(c::AbstractDigraph) -> Bool

Test whether the given AbstractDigraph `c` has children.
"""
@inline haschildren(c::AbstractDigraph) = !isempty(c.children)


"""
    isparent(p::AbstractDigraph, c::AbstractDigraph) -> Bool

Test whether `p` is the parent of `c`.
"""
@inline isparent(p::AbstractDigraph, c::AbstractDigraph) = p===c.parent
@inline isparent(::Nothing, c::AbstractDigraph) = false


"""
    setparent!(c::T, p::T) where {T<:AbstractDigraph}

Set `p` as the parent of `c`.
"""
function setparent!(child::T, parent::T) where {T<:AbstractDigraph}
    hasparent(child) && error("unable to setparent! of non-orphan item")
    push!(parent.children, child)
    child.parent = parent
    child
end


"""
    popparent!(c::AbstractDigraph}

Remove the parent from `c` and make `c` orphan.
"""
function popparent!(child::AbstractDigraph)
    if hasparent(child)
        parent = child.parent
        i = findfirst(x->x===child, parent.children)
        if i !== nothing
            deleteat!(parent.children, i)
            child.parent = nothing
        end
    end
    child
end

"""
    popchild!(p::AbstractDigraph} -> child

Remove a child from `p` and return it. The returned element is orphan.
Return `nothing` if the item has no children.
"""
popchild!(p::AbstractDigraph) = isempty(p.children) ? nothing : popparent!(p.children[1])


"""
    adjacency(p::AbstractDigraph)

Build and return the adjacency list for all nodes reachable from `p`.
"""
function adjacency(p::AbstractDigraph)
    while hasparent(p)
        p = p.parent
    end
    
    adjlist = Dict{Int, Vector{Int}}()
    bfs(p) do c
        hasparent(c) && push!(get!(adjlist, c.parent.index, []), c.index)
    end
    adjlist
end


# function node(io::IO, c::AbstractDigraph)
#     print(io, "$c(<- $(repr(c.parent))) with $(length(c.children)) children")
# end
# node(c::AbstractDigraph) = node(stdout, c)

"""
    function bfs(f::Function, c::T) where {T <: AbstractDigraph} -> Nothing

Perform a breadth-first search starting from `c` and call function
`f` on each visited node.
"""
function bfs(f::Function, c::T) where {T <: AbstractDigraph}
    queue = Vector{T}()
    push!(queue, c)
    while !isempty(queue)
        p = popfirst!(queue)
        f(p)
        append!(queue, p.children)
    end
    nothing
end
