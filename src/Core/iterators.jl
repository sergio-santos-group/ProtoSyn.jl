export eachatom, eachresidue

#region Iterators --------------------------------------------------------------
const _ByAtom = Val{1}
const _ByResidue = Val{2}

struct ItemIterator{T <: AbstractContainer,B}
    target::T
    size::Tuple{Vararg{Int}}
end

Base.length(iter::ItemIterator{T,B}) where {T,B} = iter.size[end]
Base.size(iter::ItemIterator{T,B}) where {T,B} = iter.size

Base.show(io::IO,iter::ItemIterator{T,_ByAtom}) where T =
    println(io, "ItemIterator{$(nameof(T)), _ByAtom} with size $(iter.size)")

Base.show(io::IO,iter::ItemIterator{T,_ByResidue}) where T =
    println(io, "ItemIterator{$(nameof(T)), _ByResidue} with size $(iter.size)")

# Base.show(io::IO,iter::ItemIterator{T,_ByResidue}) where T = begin
#     println(io, "ItemIterator{$(nameof(T)), _ByResidue} with size $(iter.size)")
# end

#region AtomIterator -----------------------------------------------------------

export eachatom
eachatom(c::T) where {T<:AbstractContainer} = ItemIterator{T,_ByAtom}(c, size(c))


Base.iterate(iter::ItemIterator{Topology,_ByAtom}, (s,r,a)=(1,1,1)) = begin
    # if s > length(t.segments)
    #     return nothing
    # elseif r > length(t.segments[s].residues)
    #     return iterate(iter, (s+1,1,1))
    # elseif a > length(t.segments[s].residues[r].atoms)
    #     return iterate(iter, (s,r+1,1))
    # end
    t = iter.target
    s > length(t) && return nothing
    r > length(t.items[s]) && return iterate(iter, (s+1,1,1))
    a > length(t.items[s].items[r]) && return iterate(iter, (s,r+1,1))
    (t.items[s].items[r].items[a], (s,r,a+1))
end

# Base.iterate(iter::ItemIterator{Segment,_ByAtom}, (r,a)=(1,1)) = begin
#     s = iter.target
#     if r > length(s.residues)
#         return nothing
#     elseif a > length(s.residues[r].atoms)
#         return iterate(iter, (r+1,1))
#     end
#     (s.residues[r].atoms[a], (r,a+1))
# end

Base.iterate(iter::ItemIterator{Residue,_ByAtom}, (a,)=(1,)) = begin
    r = iter.target
    if a > length(r)
        return nothing
    end
    (r.items[a], (a+1,))
end

#endregion AtomIterator


#region ResidueIterator --------------------------------------------------------
export eachresidue

eachresidue(c::T) where {T<:AbstractContainer} = 
    ItemIterator{T,_ByResidue}(c, size(c)[1:end-1])


Base.iterate(iter::ItemIterator{Topology,_ByResidue}, (s,r)=(1,1)) = begin
    t = iter.target
    if s > length(t)
        return nothing
    elseif r > length(t.items[s])
        return iterate(iter, (s+1,1))
    end
    (t.items[s].items[r], (s,r+1))
end
#region ResidueIterator


#endregion Iterators