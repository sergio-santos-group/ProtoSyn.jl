# This file should contain overloads of Julia's Base module.

#region show -------------------------------------------------------------------

name(io::IO, item::AbstractContainer) = begin
    if hascontainer(item)
        name(io, item.container)
    end
    print(io, "/$(item.name):$(item.id)")
end

Base.show(io::IO, item::AbstractContainer) = begin
    print(io, "$(nameof(typeof(item))){")
    name(io, item)
    print(io, "}")
end

#endregion show

#region push! ------------------------------------------------------------------

_push!(container::AbstractContainer{T}, item::T) where {T<:AbstractContainer} = begin
    if !in(item, container)
        push!(container.items, item)

        # The next line makes it so that a fragment can only be appended to 1
        # pose, 1 time only.
        item.container = container

        # Update the container size to reflect the addition
        container.size += 1
    end
    container
end

Base.push!(container::AbstractContainer{T}, item::T) where {T<:AbstractContainer} = begin
    _push!(container, item)
end

Base.push!(res::Residue, atm::Atom) = begin
    res.itemsbyname[atm.name] = atm
    if !in(atm, res)
        _push!(res, atm)
    end
    res
end

#endregion push!

#region insert! ----------------------------------------------------------------

Base.insert!(container::AbstractContainer{T}, index::Integer, item::T) where {T <: AbstractContainer} = begin
    insert!(container.items, index, item)
    item.container = container
    container.size += 1
    container
end

Base.insert!(container::AbstractContainer{T}, index::Integer, items::Vector{T1}) where {T <: AbstractContainer, T1 <: AbstractContainer} = begin
    for item in Iterators.reverse(items)
        insert!(container.items, index, item)
        item.container = container
        container.size += 1
    end
    container
end

Base.insert!(container::Residue, index::Integer, item::Atom) = begin
    insert!(container.items, index, item)
    item.container = container
    container.size += 1
    container.itemsbyname[item.name] = item
    container
end

#endregion insert!

#region get --------------------------------------------------------------------

Base.get(res::Residue, name::AbstractString, default = nothing) = begin
    get(res.itemsbyname, name, default)
end

Base.getindex(r::Residue, n::AbstractString) = get(r, n)

Base.getindex(c::AbstractContainer, i::Int) = c.items[i]

Base.getindex(c::AbstractContainer, i::Union{Int, AbstractString}...) = begin
    getindex(c.items[i[1]], i[2:end]...)
end

#endregion get

#region in ---------------------------------------------------------------------
@inline Base.in(item::T, c::AbstractContainer{T}) where {T<:AbstractContainer} = begin
    item.container===c
end

Base.in(name::String, r::Residue) = haskey(r.itemsbyname, name)

#endregion in

#region size -------------------------------------------------------------------

@inline Base.size(r::Residue) = (r.size,)

@inline Base.size(s::Segment) = (s.size, mapreduce(r->r.size, +, s.items; init=0))

@inline Base.size(t::Topology) = begin
    nres = natm = 0
    for seg in t.items
        nr,na = size(seg)
        nres += nr
        natm += na
    end
    (t.size, nres, natm)
end

#endregion size

#region copy -------------------------------------------------------------------

@inline Base.copy(a0::Atom) = begin
    return Atom(a0.name, a0.id, a0.index, a0.symbol)
end


@inline Base.copy(r0::Residue) = begin
    byatom = eachatom(r0)
    
    # Create residue (Automatically updates residue.index)
    r1 = Residue(r0.name, r0.id)

    # Populate residue with atoms
    # (Automatically updates residue.size, residue.items, residue.itemsbyname)
    # (Automatically updates atom.container)
    for at0 in byatom # Looping "old" residue
        push!(r1, copy(at0))
    end
    
    for at0 in byatom # Looping "old" residue
        # Build intra-residue bonds
        at1 = get(r1, at0.name) # Uses residue.itemsbyname
        for other0 in at0.bonds
            in(other0, r0) && push!(at1.bonds, get(r1, other0.name))
        end
        
        # Build intra-residue graph (set parents and children)
        if hasparent(at0)
            in(at0.parent, r0) && setparent!(at1, r1[at0.parent.name])
        end
    end

    return r1
end


Base.copy(s0::Segment) = begin
    s1 = Segment(s0.name, s0.id)
    old2new = IdDict{AbstractContainer, AbstractContainer}()
    root = origin(s0)

    # copy residues and atoms
    for r0 in eachresidue(s0)
        r1 = Residue!(s1, r0.name, r0.id)
        old2new[r0] = r1
        for at0 in eachatom(r0)
            old2new[at0] = Atom!(r1, at0.name, at0.id, at0.index, at0.symbol)
            r1.items[end].ascendents = at0.ascendents
        end
    end

    # deal with atom graph and bonds
    for at0 in eachatom(s0)
        at1 = old2new[at0]
        if hasparent(at0) && !isparent(root, at0)
            setparent!(at1, old2new[at0.parent])
        end
        for other0 in at0.bonds
            push!(at1.bonds, old2new[other0])
        end
    end

    # deal with residue graph
    for r0 in eachresidue(s0)
        hasparent(r0) && r0.parent != root.container && setparent!(old2new[r0], old2new[r0.parent])
    end

    return s1
end

Base.copy(t0::Topology) = begin
    t1 = Topology(t0.name, t0.id)
    for s0 in eachsegment(t0)
        s1 = copy(s0)
        push!(t1, s1)
        setparent!(s1[1], t1.root)
        setparent!(s1[1][1], origin(t1))
    end

    return t1
end

Base.copy(s::State{T}) where T = begin
    ns = State(T, s.size)
    # Updating item.t also updates the parent.x matrix
    for (index, atomstate) in enumerate(s)
        ns[index].t = copy(atomstate.t)
        ns[index].r = copy(atomstate.r)
        ns[index].b = atomstate.b
        ns[index].θ = atomstate.θ
        ns[index].ϕ = atomstate.ϕ
        ns[index].Δϕ = atomstate.Δϕ
        ns[index].changed = atomstate.changed
    end
    ns.id = s.id
    ns.i2c = s.i2c
    ns.c2i = s.c2i
    ns.f = copy(s.f)
    ns.e = copy(s.e)
    ns
end

Base.copy(p::Pose) = Pose(copy(p.graph), copy(p.state))

#endregion copy

#region delete! ----------------------------------------------------------------

Base.delete!(container::AbstractContainer{T}, item::T) where {T<:AbstractContainer} = begin
    if in(item, container)
        i = findfirst(x -> x === item, container.items)
        if i !== nothing
            deleteat!(container.items, i)
            item.container = nothing
            container.size -= 1
        end
    end
    container
end

#endregion delete!

#region Others ----------------------------------------------------------------

@inline Base.isempty(c::AbstractContainer) = c.size==0

@inline Base.length(c::AbstractContainer) = c.size

@inline Base.lastindex(c::AbstractContainer) = c.size

@inline Base.firstindex(c::AbstractContainer) = 1

@inline Base.haskey(r::Residue, k::AbstractString) = haskey(r.itemsbyname, k)

# @inline Base.findfirst(x::T, c::Vector{T}) where T = findfirst(i -> i === x, c) # ???

@inline Base.findfirst(x::T, c::AbstractContainer{T}) where T = findfirst(x, c.items)

#endregion Others

#region pop! -------------------------------------------------------------------

Base.pop!(pose::Pose{Topology}, seg::Segment) = begin

    if seg.container !== pose.graph
        error("given residue does not belong to the provided topology")
    end

    # remove node states from parent state and create
    # a new state for this segment
    st = splice!(pose.state, seg[1][1].index:seg[end][end].index)
    
    pop!(pose.graph.items, seg)

    reindex(pose.graph)
    
    # new common ID
    seg.id = st.id = genid()
    Pose(seg, st)
end

Base.pop!(top::Topology, seg::Segment) = begin
    deleteat!(top.items, findall(x -> x == seg, top.items))
    top.size -= 1
end

#endregion pop!

#region detach -----------------------------------------------------------------

Base.detach(s::Segment) = begin
    root = origin(s)
    for at in root.children
        isparent(root, at) && popparent!(at)
    end
    hascontainer(s) && delete!(s.container, s)
end

#endregion detach