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

    return container
end

"""
    push!(container::AbstractContainer{T}, item::T) where {T <: AbstractContainer}

Add an `AbstractContainer` `item` to the `AbstractContainer` `container`,
updating the `container` size and setting the correct `item.container`. Return
the altered `container`.

    push!(residue::Residue, atom::Atom)

In the specific case of adding an [`Atom`](@ref) `atom` to a [`Residue`](@ref)
`residue`, also add the `atom.name` to the `residue.itemsbyname` dictionary for
correct indexation by name. Return the altered [`Residue`](@ref) `residue`.

    push!(container::AbstractContainer{T}, items::Vector{T}) where {T <: AbstractContainer}

Add a vector of `AbstractContainer` `items` in the `AbstractContainer`
`container`, updating the `container` size and setting the correct
`item.container` for each item in the `items`. Note that this method keeps the
order of `items`. Return the altered `container`.

*This function is a Base module overload.*

!!! ukw "Note:"
    This function does not set any [Bonds](@ref) or
    [Parenthood relationships](@ref) to other items in the same `container`.
    This function does not set a complementary [State](@ref state-types).

# See also
[`insert!`](@ref) [`delete!`](@ref)

# Examples
```
julia> push!(pose.graph[1][1], Atom("CA", -1, -1, "C"))
Residue{/UNK:1/UNK:1/SER:1}
```
"""
Base.push!(container::AbstractContainer{T}, item::T) where {T <: AbstractContainer} = begin
    return _push!(container, item)
end

Base.push!(container::AbstractContainer{T}, items::Vector{T}) where {T <: AbstractContainer} = begin
    for item in items
        _push!(container, item)
    end

    return container
end

Base.push!(residue::Residue, atom::Atom) = begin
    if atom.name in keys(residue.itemsbyname)
        @warn "Atom named $(atom.name) already found in residue $residue. Please check for the existance of duplicates."
    end
    residue.itemsbyname[atom.name] = atom
    if !in(atom, residue)
        _push!(residue, atom)
    end

    return residue
end

#endregion push!

#region insert! ----------------------------------------------------------------

"""
    insert!(container::AbstractContainer{T}, index::Integer, item::T) where {T <: AbstractContainer}

Insert an `AbstractContainer` `item` in the `AbstractContainer` `container` at
the given `index`, updating the `container` size and setting the correct
`item.container`. Return the altered `container`.

    insert!(container::Residue, index::Integer, item::Atom)

In the specific case of inserting an [`Atom`](@ref) `atom` in a
[`Residue`](@ref) `residue`, also add the `atom.name` to the
`residue.itemsbyname` dictionary for correct indexation by name. Return the
altered [`Residue`](@ref) `residue`.

    insert!(container::AbstractContainer{T}, index::Integer, items::Vector{T}) where {T <: AbstractContainer}

Insert a vector of `AbstractContainer` `items` in the `AbstractContainer`
`container` at the given `index`, updating the `container` size and setting the
correct `item.container` for each item in the `items`. Note that this method
keeps the order of `items`. Return the altered `container`.

*This function is a Base module overload.*

!!! ukw "Note:"
    This function does not set any [Bonds](@ref) or
    [Parenthood relationships](@ref) to other items in the same `container`.
    This function does not set a complementary [State](@ref state-types).

# See also
[`push!`](@ref) [`delete!`](@ref)

# Examples
```jldoctest
julia> insert!(pose.graph[1][1], 3, Atom("CA", -1, -1, "C"))
Residue{/UNK:1/UNK:1/SER:1}
```
"""
Base.insert!(container::ProtoSyn.AbstractContainer{T}, index::Integer, item::T) where {T <: AbstractContainer} = begin
    insert!(container.items, index, item)
    item.container = container
    container.size += 1

    return container
end

Base.insert!(container::ProtoSyn.AbstractContainer{T}, index::Integer, items::Vector{T}) where {T <: AbstractContainer} = begin
    for item in Iterators.reverse(items)
        insert!(container.items, index, item)
        item.container = container
        container.size += 1
    end

    return container
end

Base.insert!(container::Residue, index::Integer, item::Atom) = begin
    insert!(container.items, index, item)
    item.container = container
    container.size += 1
    container.itemsbyname[item.name] = item

    return container
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

"""
    copy(container::AbstractContainer)

Return a copy of the given `AbstractContainer` `container`. Note that for
[`Residue`](@ref) instances only the intra-residue [Bonds](@ref) and
[Parenthood relationships](@ref) can be copied. The same logic applies to the
intra-segment [Graph](@ref graph-types).

# Example
```jldoctest
julia> copy(pose.graph[1, 1, 1])
Atom{/N:1}

julia> copy(pose.graph[1, 1])
Residue{/SER:1}
```
"""
@inline Base.copy(a0::Atom) = begin
    return Atom(a0.name, a0.id, a0.index, a0.symbol)
end


@inline Base.copy(r0::Residue) = begin
    byatom = eachatom(r0)
    
    # Create residue (Automatically updates residue.index)
    r1 = Residue(r0.name.content, r0.id)

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
    s1      = Segment(s0.name, s0.id)
    s1.code = s0.code
    old2new = IdDict{AbstractContainer, AbstractContainer}()
    root    = ProtoSyn.root(s0)

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
    # for at0 in ProtoSyn.travel_graph(s0[1, 1]) # This option doesn't set all relationships
        at1 = old2new[at0]
        if hasparent(at0) && !isparent(root, at0) && at0.parent in keys(old2new)
            setparent!(at1, old2new[at0.parent])
        end
        for other0 in at0.bonds
            other0 in keys(old2new) && push!(at1.bonds, old2new[other0])
        end
    end

    # deal with residue graph
    for r0 in eachresidue(s0)
        if root === nothing # Dealing with fragment
            hasparent(r0) && setparent!(old2new[r0], old2new[r0.parent])
        else
            hasparent(r0) && r0.parent != root.container && setparent!(old2new[r0], old2new[r0.parent])
        end
    end

    return s1
end

Base.copy(t0::Topology) = begin
    t1 = Topology(t0.name, genid())
    for s0 in eachsegment(t0)
        s1 = copy(s0)
        push!(t1, s1)
        # setparent!(s1[1], t1.root)
        # setparent!(s1[1][1], root(t1))
    end

    root0 = ProtoSyn.root(t0)
    root1 = ProtoSyn.root(t1)
    for (atom0, atom1) in zip(eachatom(t0), eachatom(t1))
        if atom0.parent == root0
            setparent!(atom1, root1)
        end
    end

    for (res0, res1) in zip(eachresidue(t0), eachresidue(t1))
        if res0.parent == root0.container
            setparent!(res1, root1.container)
        end
    end

    return t1
end

Base.copy(s::State{T}) where T = begin
    ns = State(T, s.size)
    # Updating item.t also updates the parent.x matrix
    
    # Also copies the current root coordinates, if they are present
    if any((as) -> as.index < 0, s.items)
        for (index, root_atomstate) in enumerate(s.items[1:3])
            ns[index - 3].t = copy(root_atomstate.t)
            ns[index - 3].r = copy(root_atomstate.r)
            ns[index - 3].b = root_atomstate.b
            ns[index - 3].θ = root_atomstate.θ
            ns[index - 3].ϕ = root_atomstate.ϕ
            ns[index - 3].Δϕ = root_atomstate.Δϕ
            ns[index - 3].changed = root_atomstate.changed
            ns[index - 3].δ = root_atomstate.δ
        end
    end

    for (index, atomstate) in enumerate(s)
        ns[index].t = copy(atomstate.t)
        ns[index].r = copy(atomstate.r)
        ns[index].b = atomstate.b
        ns[index].θ = atomstate.θ
        ns[index].ϕ = atomstate.ϕ
        ns[index].Δϕ = atomstate.Δϕ
        ns[index].changed = atomstate.changed
        ns[index].δ = atomstate.δ
    end
    ns.id = genid()
    ns.i2c = s.i2c
    ns.c2i = s.c2i
    ns.f = copy(s.f)
    ns.e = copy(s.e)
    ns
end


"""
    copy(pose::Pose)

Return a copied [Pose](@ref pose-types) of the provided `pose`. The resulting [Pose](@ref pose-types)
will have different `:id` fields for the [Graph](@ref state-types) [`Topology`](@ref) and
[State](@ref state-types).

*This function is a Base module overload.*

# Examples
```
julia> copy(pose)
Pose{Topology}(Topology{/UNK:9547}, State{Float64}:
 Size: 343
 i2c: false | c2i: false
 Energy: Dict(:Total => Inf)
)
```
"""
Base.copy(pose::Pose) = begin
    graph = copy(pose.graph)
    state = copy(pose.state)
    graph.id = state.id = genid()
    return Pose(graph, state)
end

#endregion copy

#region delete! ----------------------------------------------------------------

"""
    delete!(container::AbstractContainer{T}, item::T) where {T <: AbstractContainer}

Delete the given `AbstractContainer` `item` from the `AbstractContainer`
`container`, if found, while updating the `container.size` and `item.container`
fields. In the specific case of deleting an [`Atom`](@ref) instance from a
[`Residue`](@ref), update `container.itemsbyname` as well. Return the altered
`container`. If the given `AbstractContainer` `item` is not found in the
`AbstractContainer` `container`, return `nothing`.

# See also
[`push!`](@ref) [`insert!`](@ref)

# Examples
```jldoctest
julia> delete!(pose.graph[1][1], pose.graph[1][1]["CA"])
Residue{/UNK:1/UNK:1/SER:1}
```
"""
Base.delete!(container::AbstractContainer{T}, item::T) where {T <: AbstractContainer} = begin
    if in(item, container)
        i = findfirst(x -> x === item, container.items)
        if i !== nothing
            if hasproperty(container, :itemsbyname)
                delete!(container.itemsbyname, item.name)
            end
            deleteat!(container.items, i)
            item.container = nothing
            container.size -= 1
        end
    end

    return container
end

Base.delete!(container::AbstractContainer{T}, item::Nothing) where {T <: AbstractContainer}= begin
    @warn "Tried to delete an instance of $T which was not found in the given container $container."
    return item
end

#endregion delete!

#region Others ----------------------------------------------------------------

@inline Base.isempty(c::AbstractContainer) = c.size==0

@inline Base.length(c::AbstractContainer) = c.size

@inline Base.lastindex(c::AbstractContainer) = c.size

@inline Base.firstindex(c::AbstractContainer) = 1

@inline Base.haskey(r::Residue, k::AbstractString) = haskey(r.itemsbyname, k)

@inline Base.findfirst(x::T, c::Vector{T}) where {T <: AbstractContainer} = findfirst(i -> i === x, c)

@inline Base.findfirst(x::T, c::AbstractContainer{T}) where T = findfirst(x, c.items)

#endregion Others

#region detach -----------------------------------------------------------------

"""

    detach(segment::Segment)

Detach and return the given [`Segment`](@ref) from it's container
[Graph](@ref state-types), by:
* Detaching any [`Atom`](@ref) and [`Residue`](@ref) instance from the [Graph](@ref state-types)'s Root (by popping parenthood relationships), if said instances belong to the given [`Segment`](@ref) instance.
* Deleting this [`Segment`](@ref) from its container [`Topology`](@ref). 

*This function is a Base module overload.*

!!! ukw "Note:"
    This function does not alter the [State](@ref state-types) of the [Pose](@ref pose-types)
    containing the provided [`Segment`](@ref). 

# Examples
```jldoctest
julia> detach(pose.graph[1])
Segment{/UNK:1}
```
"""
Base.detach(segment::Segment) = begin
    # * This function should remove the parenthood to the origin on all atoms
    # * and residues of this segment that are connected to the root.
    root = ProtoSyn.root(segment)
    for atom in root.children
        !(atom.container.container == segment) && continue
        isparent(root, atom) && popparent!(atom)
        isparent(root.container, atom.container) && popparent!(atom.container)
    end
    hascontainer(segment) && delete!(segment.container, segment)
    return segment
end

#endregion detach