export link!, unlink!
export hasparent, haschildren


@inline hasparent(c::AbstractContainer) = c.parent !== nothing
@inline haschildren(c::AbstractContainer) = !isempty(c.children)


@doc """
set parent of child
"""
function link!(parent::T, child::T) where {T<:AbstractContainer}
# function setparent!(parent::T, child::T) where {T<:AbstractContainer}
    hasparent(child) && error("unable to link! of non-orphan item (perhaps it if a child of root?)")

    push!(parent.children, child)
    child.parent = parent
    parent
end


#Base.delete!(parent::T, child::T) where {T<:AbstractContainer} = begin
function unlink!(parent::T, child::T) where {T<:AbstractContainer}
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
    print(io, "$c(<- $(repr(c.parent))) with $(length(c.children)) children")
end
node(c::AbstractContainer) = node(stdout, c)
