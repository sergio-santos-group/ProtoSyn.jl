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

Base.insert!(container::AbstractContainer{T}, index::Integer, items::Vector{T}) where {T <: AbstractContainer} = begin
    for item in Iterators.reverse(iterator(T)(items))
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

@inline Base.copy(a::Atom) = begin
    Atom(a.name, a.id, a.index, a.symbol)
end

@inline Base.copy(r::Residue) = begin
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
        if hasparent(at)
            in(at, r) && setparent!(at1, r1[at.parent.name])
        end
    end
    r1
end

Base.copy(s::Segment) = begin
    old2new = IdDict{AbstractContainer,AbstractContainer}()
    s1 = Segment(s.name, s.id)
    root = origin(s)

    # copy residues and atoms
    for r in eachresidue(s)
        r1 = Residue!(s1, r.name, r.id)
        old2new[r] = r1
        for at in eachatom(r)
            old2new[at] = Atom!(r1, at.name, at.id, at.index, at.symbol)
        end
    end

    # deal with atom graph and bonds
    for at in eachatom(s)
        at1 = old2new[at]
        hasparent(at) && !isparent(root, at) && setparent!(at1, old2new[at.parent])
        for other in at.bonds
            push!(at1.bonds, old2new[other])
        end
    end

    # deal with residue graph
    for r in eachresidue(s)
        hasparent(r) && setparent!(old2new[r], old2new[r.parent])
    end

    s1
end

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

@inline Base.findfirst(x::T, c::Vector{T}) where T = findfirst(i -> i === x, c)

@inline Base.findfirst(x::T, c::AbstractContainer{T}) where T = findfirst(x, c.items)

#endregion Others