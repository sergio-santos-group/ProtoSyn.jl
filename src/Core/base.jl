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
        item.container = container
        container.size += 1
    end
    container
end

Base.push!(container::AbstractContainer{T}, item::T) where {T<:AbstractContainer} = begin
     _push!(container, item)
end

Base.push!(res::Residue, atm::Atom) = begin
    if !in(atm,res)
        res.itemsbyname[atm.name] = atm
        _push!(res, atm)
    end
    res
end

Base.insert!(container::AbstractContainer{T}, index::Integer, item::T) where {T<:AbstractContainer} = begin
    insert!(container.items, index, item)
    item.container = container
    container.size += 1
    container
end

#endregion push!


#region get --------------------------------------------------------------------


Base.get(res::Residue, name::String, default=nothing) = begin
    get(res.itemsbyname, name, default)
end

Base.getindex(r::Residue, n::AbstractString) = get(r,n)
Base.getindex(c::AbstractContainer, i::Int) = c.items[i]
Base.getindex(c::AbstractContainer, i::Union{Int,AbstractString}...) = begin
    getindex(c.items[i[1]], i[2:end]...)
end


#endregion get


@inline Base.in(item::T, c::AbstractContainer{T}) where {T<:AbstractContainer} = begin
    item.container===c
end


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

@inline Base.isempty(c::AbstractContainer) = c.size==0
@inline Base.length(c::AbstractContainer) = c.size

@inline Base.lastindex(c::AbstractContainer) = c.size

@inline Base.firstindex(c::AbstractContainer) = 1



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


Base.delete!(container::AbstractContainer{T}, item::T) where {T<:AbstractContainer} = begin
    if in(item, container)
        i = findfirst(x->x===item, container.items)
        if i !== nothing
            deleteat!(container.items, i)
            item.container = nothing
            container.size -= 1
        end
    end
    container
end





export hascontainer
@inline hascontainer(c::AbstractContainer) = c.container !== nothing

# QUESTION
# Whatever we do, origin(something) is always the ROOT? In what case for this
# function return 'nothing'?

@inline origin(t::Topology) = get(t.root, "OO")
@inline origin(c::AbstractContainer) = hascontainer(c) ? origin(c.container) : nothing

@inline genid() = Int(rand(UInt16))

export reindex
# @inline reindex(c::AbstractContainer) = begin
@inline reindex(t::Topology) = begin
    # index = 0
    aid = rid = sid = 0
    for seg in t.items
        seg.id = (sid +=1 )
        for res in seg.items
            res.id = (rid +=1 )
            for atm in res.items
                atm.index = (aid +=1 )
            end
        end
    end
    # update ascendents (not possible before because
    # of possible problems with index assignment)
    for atm in eachatom(t)
        println(atm, atm === nothing)
        atm.ascendents = ascendents(atm, 4)
    end

    t
end

reindex(s::Segment) = begin
    aid = rid = 0
    for res in s.items
        res.id = (rid += 1)
        for atm in res.items
            atm.index = (aid +=1 )
        end
    end
    s
end

Base.haskey(r::Residue, k::AbstractString) = haskey(r.itemsbyname, k)