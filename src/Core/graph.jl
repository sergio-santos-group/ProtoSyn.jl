export setparent!
function setparent!(parent::T, child::T) where {T<:AbstractContainer}
    if parent.name=="O"
        println("aqui", parent, child)
    end
    push!(parent.children, child)
    child.parent = parent
    parent
end

export hasparent
@inline hasparent(c::AbstractContainer) = c.parent !== nothing




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