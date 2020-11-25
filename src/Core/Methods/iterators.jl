export eachatom, eachresidue

#region Iterators --------------------------------------------------------------
# const _ByAtom = Val{1}
# const _ByResidue = Val{2}
struct _ByAtom end
struct _ByResidue end
struct _BySegment end

struct ItemIterator{T <: AbstractContainer, B}
    target::T
    size::Tuple{Vararg{Int}}
end

Base.size(iter::ItemIterator{T,B}) where {T,B} = iter.size
Base.length(iter::ItemIterator{T,B}) where {T,B} = iter.size[end]

Base.show(io::IO,iter::ItemIterator{T,B}) where {T,B} =
    println(io, "ItemIterator{$(nameof(T)), $(nameof(B))} with size $(iter.size)")

#region AtomIterator -----------------------------------------------------------

export eachatom
eachatom(c::T) where {T<:AbstractContainer} = ItemIterator{T, _ByAtom}(c, size(c))


Base.iterate(iter::ItemIterator{Topology, _ByAtom}, (s, r, a)=(1, 1, 1)) = begin
    t = iter.target
    s > length(t) && return nothing
    r > length(t.items[s]) && return iterate(iter, (s+1,1,1))
    a > length(t.items[s].items[r]) && return iterate(iter, (s,r+1,1))
    (t.items[s].items[r].items[a], (s,r,a+1))
end

Base.iterate(iter::ItemIterator{Segment, _ByAtom}, (r,a)=(1,1)) = begin
    s = iter.target
    if r > length(s.items)
        return nothing
    elseif a > length(s.items[r].items)
        return iterate(iter, (r+1,1))
    end
    (s.items[r].items[a], (r,a+1))
end

Base.iterate(iter::ItemIterator{Residue, _ByAtom}, (a,)=(1,)) = begin
    r = iter.target
    if a > length(r)
        return nothing
    end
    (r.items[a], (a+1,))
end

#endregion AtomIterator


#region ResidueIterator --------------------------------------------------------
export eachresidue

eachresidue(c::T) where {T <: AbstractContainer} = 
    ItemIterator{T, _ByResidue}(c, size(c)[1:end-1])


Base.iterate(iter::ItemIterator{Topology, _ByResidue}, (s, r)=(1, 1)) = begin
    t = iter.target
    if s > length(t)
        return nothing
    elseif r > length(t.items[s])
        return iterate(iter, (s+1, 1))
    end
    (t.items[s].items[r], (s, r+1))
end

Base.iterate(iter::ItemIterator{Segment, _ByResidue}, (r,)=(1,)) = begin
    s = iter.target
    if r > length(s)
        return nothing
    end
    (s.items[r], (r+1,))
end

#region ResidueIterator

# QUESTION: This works?
export eachsegment
eachsegment(t::Topology) = ItemIterator{Topology, _BySegment}(t, (length(t),))
Base.iterate(iter::ItemIterator{Topology, _BySegment}, (s,)=(1,)) = begin
    t = iter.target
    if s > length(t)
        return nothing
    end
    (t.items[s], (s+1,))
end

export counter
counter(::Type{Segment})  = count_segments
counter(::Type{Residue})  = count_residues
counter(::Type{Atom})     = count_atoms

export iterator
iterator(::Type{Segment}) = eachsegment
iterator(::Type{Residue}) = eachresidue
iterator(::Type{Atom})    = eachatom



#endregion Iterators