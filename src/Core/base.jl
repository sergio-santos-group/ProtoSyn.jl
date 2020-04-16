#region show -------------------------------------------------------------------

name(io::IO, item::AbstractContainer) = begin
    if isdefined(item, :container) && item.container !== nothing
        name(io, item.container)
    end
    print(io, "/$(item.name):$(item.id)")
end

Base.show(io::IO, item::AbstractContainer) = begin
    print(io, "$(nameof(typeof(item))){")
    name(io, item)
    print(io, "}")
end

# Base.show(io::IO, item::T) where {T<:AbstractContainer}= begin
#     print(io, "$(nameof(T)){")
#     name(io, item)
#     print(io, "}")
# end

# Base.show(io::IO, node::GraphNode{T}) where T = begin
#     parent = node.parent===nothing ? "nothing" : node.parent.item
#     println("GraphNode{$(node.item), $parent} with $(length(node.children)) children")
# end

#endregion show

@inline Base.in(item::T, c::AbstractContainer{T}) where {T<:AbstractContainer} = begin
    item.container===c
end

#region push! ------------------------------------------------------------------

_push!(container::AbstractContainer{T}, item::T) where {T<:AbstractContainer} = begin
    if !in(item, container)
        push!(container.items, item)
        item.container = container
        container.size += 1
    end
    container
end

# Base.push!(container::AbstractContainer{T}, item::T) where {T<:AbstractContainer} = begin
Base.push!(container::AbstractContainer{T}, item::T) where T = _push!(container, item)

Base.push!(res::Residue, atm::Atom) = begin
    if !in(atm,res)
        res.itemsbyname[atm.name] = atm
        _push!(res, atm)
    end
    res
end


# Base.push!(parent::GraphNode{T}, node::GraphNode{T}) where T = begin
#     push!(parent.children, node)
#     node.parent = parent
#     parent
# end



#endregion push!







#region get --------------------------------------------------------------------
Base.getindex(r::Residue, n::AbstractString) = get(r,n)
Base.get(res::Residue, name::String, default=nothing) = 
    get(res.itemsbyname, name, default)

Base.getindex(c::AbstractContainer, i::Int) = c.items[i]
Base.getindex(c::AbstractContainer, i::Int...) = getindex(c.items[i[1]],i[2:end]...)

#Base.get(res::Residue, name::String) = get(res, name, nothing)

#endregion get


# export  isorphan
# @inline isorphan(node::GraphNode) = begin
#     node.parent===nothing && isempty(node.children)
# end

# @inline isorphan(c::AbstractContainer) = c.container===nothing


@inline Base.in(item::T, container::AbstractContainer{T}) where T = begin
    item.container===container
end


# export isorphan
# @inline isorphan(c::AbstractContainer) = c.container === nothing

export hascontainer
@inline hascontainer(c::AbstractContainer) = c.container !== nothing




# @inline hasparent(c::AbstractContainer) = c.container !== nothing
# @inline hasparent(n::GraphNode) = n.parent !== nothing
# # @inline hasparent(::Nothing) = false

@inline origin(t::Topology) = get(t.root, "OO")
@inline origin(c::AbstractContainer) = 
    hascontainer(c) ? origin(c.container) : nothing

@inline genid() = Int(rand(UInt16))

export reindex
@inline reindex(c::AbstractContainer) = begin
    index = 0
    for atom in eachatom(c)
        atom.index = (index += 1)
    end
    c
end

Base.size(r::Residue) = begin
    (r.size,)
end
Base.size(s::Segment) = begin
    (s.size, mapreduce(r->r.size, +, s.items; init=0))
end
Base.size(t::Topology) = begin
    nseg = nres = 0
    for seg in t.items
        ns,nr = size(seg)
        nseg += ns
        nres += nr
    end
    (t.size, nseg, nres)
end

Base.isempty(c::AbstractContainer) = c.size==0
Base.length(c::AbstractContainer) = c.size


# Base.getindex(c::AbstractContainer, i::Int) = c.items[i]

# Base.getindex(t::Topology, s::Int, r::Int) = t.items[s].items[r]
# Base.getindex(t::Topology, s::Int, r::Int, a::Int) = t.items[s].items[r].items[a]
# # Base.getindex(s::Topology, s::Int, r::Int) = t.items[s].items[r]

Base.lastindex(c::AbstractContainer) = c.size
Base.firstindex(c::AbstractContainer) = 1



Base.copy(a::Atom) = begin
    Atom(a.name, a.id, a.index, a.symbol)
end

Base.copy(r::Residue) = begin
    byatom = eachatom(r)
    
    # create residue and copy atoms
    r1 = Residue(r.name, r.id)
    for at in byatom
        push!(r1, copy(at))
    end
    
    for at in byatom
        # add intra-residue bonds
        at1 = get(r1, at.name)
        for other in at.bonds
            in(other,r) && push!(at1.bonds, get(r1, other.name))
        end
        
        # build atom graph
        # if (pnode = parent(at)) !== nothing
        #     in(pnode.item, r) && push!(get(r1, pnode.item.name).node, at1.node)
        # end
        if hasparent(at)
            # in(at, r) && setparent!(get(r1, at.parent.name), at1)
            in(at, r) && setparent!(r1[at.parent.name], at1)
        end
    end
    r1
end



