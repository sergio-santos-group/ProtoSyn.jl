export setparent!
export hasparent,haschildren


@inline hasparent(c::AbstractContainer) = c.parent !== nothing
@inline haschildren(c::AbstractContainer) = !isempty(c.children)


@doc """
set parent of child
"""
function setparent!(parent::T, child::T) where {T<:AbstractContainer}
    hasparent(child) && error("unable to setparent! of non-orphan item (perhaps it if a child of root?)")

    push!(parent.children, child)
    child.parent = parent
    parent
end




Base.delete!(parent::T, child::T) where {T<:AbstractContainer} = begin
    if child.parent === parent
        i = findfirst(x->x===child, parent.children)
        if i !== nothing
            deleteat!(parent.children, i)
            child.parent = nothing
        end
    end
    parent
end



function adjacency(c::AbstractContainer)
    p = c
    while hasparent(p)
        p = p.parent
    end
    
    adjlist = Dict{Int, Vector{Int}}()
    queue = AbstractContainer[]
    
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


function node(io::IO, c::AbstractContainer)
    ctype = nameof(typeof(c))
    parent = hasparent(c) ? c.parent : "nothing"
    print(io, "$ctype{$node <- $parent} with $(length(c.children)) children")
end
node(c::AbstractContainer) = node(stdout, c)
