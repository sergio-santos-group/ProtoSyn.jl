export hasparent, haschildren
export setparent!, popparent!
export popchild!

# Objects implementing the AbstractDigraph interface
# should have the following fields
#  - parent::Union{T,Nothing}
#  - children::Vector{T}
#  - visited::Bool
#  - ascendents::Union{Nothing,NTuple{N,Int}}

abstract type AbstractDigraph end


@inline hasparent(c::AbstractDigraph) = c.parent !== nothing
@inline haschildren(c::AbstractDigraph) = !isempty(c.children)


function initgraph!(c::T) where {T<:AbstractDigraph}
    c.children = Vector{T}()
    c.ascendents = nothing
    c.parent = nothing
    c.visited = false
    c
end


@doc """
set parent of child
"""
function setparent!(child::T, parent::T) where {T<:AbstractDigraph}
    hasparent(child) && error("unable to setparent! of non-orphan item")
    push!(parent.children, child)
    child.parent = parent
    child
end


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


popchild!(p::AbstractDigraph) = popparent!(pop!(p.children))


function adjacency(p::AbstractDigraph)
    while hasparent(p)
        p = p.parent
    end
    
    adjlist = Dict{Int, Vector{Int}}()
    # bfs(p) do c
    #     if !haskey(adjlist, c.parent.index)
    #         c.parent.index = []
    #     end
    #     push!(adjlist[c.parent.index], c.index)
    # end

    queue = AbstractDigraph[]
    push!(queue, p)
    while !isempty(queue)
        p = popfirst!(queue)
        lst = []
        for c in p.children
            push!(lst, c.index)
            push!(queue, c)
        end
        adjlist[p.index] = lst
    end
    adjlist
    
end


function node(io::IO, c::AbstractDigraph)
    print(io, "$c(<- $(repr(c.parent))) with $(length(c.children)) children")
end
node(c::AbstractDigraph) = node(stdout, c)


function bfs(f::Function, c::T) where {T <: AbstractDigraph}
    queue = Vector{T}()
    push!(queue, c)
    while !isempty(queue)
        p = popfirst!(queue)
        f(p)
        append!(queue, p.children)
    end
end
