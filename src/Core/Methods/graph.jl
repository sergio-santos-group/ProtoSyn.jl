# This file should contain functions that work on the system Graph, such as 
# functions that deal with parent/child relations, indexes and bonds, among
# others.

export origin
"""
    origin(container::AbstractContainer)

Return the first `Atom` in `AbstractContainer` `container` that has no parent.
The iteration follows the [`Atom`](@ref) instance `:id` field, if correctly
indexed. If no [`Atom`](@ref) instance without parent is found (i.e.: circular
structures), return `nothing`.

# See also
[`root`](@ref) [`reindex`](@ref)

"""
origin(container::AbstractContainer) = begin
    for atom in eachatom(container)
        !hasparent(atom) && return atom
    end
    return nothing
end


export root
"""
    root(container::AbstractContainer)

Return the first [`Atom`](@ref) of the Root of the [Graph](@ref state-types) that given
`AbstractContainer` `container` belongs to. If the given `AbstractContainer`
`container` is not a [`Topology`](@ref) instance and has `:container` field set
to `nothing`, return `nothing`.

    root(topology::Topology)

Return the first [`Atom`](@ref) of the Root of the given [`Topology`](@ref)
`topology` instance.

# See also
[`origin`](@ref)
"""
@inline root(topology::Topology) = get(topology.root, "OO")
@inline root(container::AbstractContainer) = begin
    hascontainer(container) ? root(container.container) : nothing
end


export hasparent
"""
    hasparent(c::AbstractContainer) -> Bool

Test whether the given AbstractContainer `c` has a parent.

# See Also
[`isparent`](@ref)
"""
@inline hasparent(c::AbstractContainer) = c.parent !== nothing


export isparent
"""
    isparent(parent::AbstractContainer, child::AbstractContainer)

Test whether `parent` is the parent of `child`.

# See Also
[`hasparent`](@ref)
"""
@inline isparent(parent::AbstractContainer, child::AbstractContainer) = begin
    parent === child.parent
end

@inline isparent(::Nothing, child::AbstractContainer) = false


export setparent!
"""
    setparent!(child::T, parent::T) where {T <: AbstractContainer}

Set `parent` as the parent of `child`, while adding `child` to
`parent.children`.

# See also
[`popparent!`](@ref)
"""
function setparent!(child::T, parent::T) where {T <: AbstractContainer}
    hasparent(child) && begin
        error("unable to setparent! of non-orphan item")
    end
    push!(parent.children, child)
    child.parent = parent
    child
end


export popparent!
"""
    popparent!(child::AbstractContainer}

Remove the parent from `child` (sets it to `nothing`) while removing `child`
from `parent.children` (only if `child` is a child of `parent`).

# See also
[`setparent!`](@ref)
"""
function popparent!(child::AbstractContainer)
    if hasparent(child)
        parent = child.parent
        i = findfirst(x -> x === child, parent.children)
        if i !== nothing
            deleteat!(parent.children, i)
            child.parent = nothing
        end
    else
        @warn "Tried to pop parent of $child but no parent was found."
    end
    child
end


export hascontainer
"""
    hascontainer(c::AbstractContainer)

Return `true` if the given `AbstractContainer.container` is not `nothing`.
"""
@inline hascontainer(c::AbstractContainer) = c.container !== nothing


export genid
"""
    genid()

Return a random `UInt16` number.
"""
@inline genid() = Int(rand(UInt16))


export reindex
"""
    reindex(topology::Topology; set_ascendents = true)
    
Re-indexes the whole [`Topology`](@ref) `topology`, setting both the `:id` and
`:index` of instances inside the `topology` to the corresponding relative index
in the `container.items` which they belong to. If `set_ascendents` is set to
`true` (is, by default), each [`Atom`](@ref) instance `:ascendents` field will
be updated to reflect the new indices.

    reindex(segment::Segment)

Re-indexes a [`Segment`](@ref) `segment`, setting both the `:id` and `:index` of
instances inside the `topology` to the corresponding relative index in the
`container.items` which they belong to.

# See also
[`ascendents`](@ref) [`reindex(::State)`](@ref)

# Examples
```jldoctest
julia> reindex(pose.graph)

julia> reindex(pose.graph[1])
```
"""
function reindex(topology::Topology; set_ascendents::Bool = true)
    aid = rid = sid = 0
    for segment in topology.items
        segment.id    = (sid += 1)
        segment.index = sid
        for residue in segment.items
            residue.id    = (rid += 1)
            residue.index = rid
            for atom in residue.items
                atom.id    = (aid += 1)
                atom.index = aid
            end
        end
    end

    # update ascendents (not possible before because
    # of possible problems with index assignment)
    if set_ascendents
        for atom in eachatom(topology)
            atom.ascendents = ascendents(atom, 4)
        end
    end

    return topology
end

function reindex(segment::Segment)
    aid = rid = 0
    for residue in segment.items
        residue.id    = (rid += 1)
        residue.index = rid
        for atom in residue.items
            atom.id    = (aid += 1)
            atom.index = aid
        end
    end
    return segment
end


export ascendents
"""
    ascedents(container::AbstractContainer, level::Int)
    
Return a `Tuple` containing the N (`level`) previous `:id` fields of the
`:parent` `AbstractContainer` instances of the given `container` (recursivelly).

# Examples
```jldoctest
julia> ascendents(atom, 4)
(1, 2, 3, 4)
```
"""
function ascendents(container::AbstractContainer, level::Int)
    if level > 1
        return (container.index, ascendents(container.parent, level - 1)...)
    else
        (container.index,)
    end
end


export unbond
"""
    unbond(pose::Pose, at1::Atom, at2::Atom)::Pose
    
Return a [Pose](@ref) instance with both given [`Atom`](@ref) instances unbonded
(removed from eachother `bonds` list; pops parenthood - if parent - and, if bond
is inter-residue, sets the downstream [`Residue`](@ref) to be the Root of the
[`Topology`](@ref). 

!!! ukw "Note:"
    This function does not maintain the  same relative position as measured from
    internal coordinates, after [`sync!`](@ref) is called. In order to keep the
    position of the downstream [`Residue`](@ref) instances, call
    [`request_c2i!`](@ref) and [`sync!`](@ref), in order to _fixate_ the
    cartesian coordinates to the new [Graph](@ref state-types) organization.

# Examples
```jldoctest
julia> unbond(pose, atom1, atom2)
```
"""
function unbond(pose::Pose, at1::Atom, at2::Atom)::Pose
    @assert (at2 in at1.bonds) & (at1 in at2.bonds) "Atoms $at1 and $at2 are not bonded and therefore cannot be unbonded."
    isparent(at1, at2) && return _unbond(pose, at1, at2)
    isparent(at2, at1) && return _unbond(pose, at2, at1)
    
    # The two atoms might be bonded but not have a parenthood relationship, in
    # which case we just remove eachother from the bond list and return the pose
    i = findfirst(at1, at2.bonds)
    i !== nothing && deleteat!(at2.bonds, i)
    j = findfirst(at2, at1.bonds)
    j !== nothing && deleteat!(at1.bonds, j)
    return pose
end

function _unbond(pose::Pose, at1::Atom, at2::Atom)::Pose
    
    i = findfirst(at1, at2.bonds)
    i !== nothing && deleteat!(at2.bonds, i)
    
    j = findfirst(at2, at1.bonds)
    j !== nothing && deleteat!(at1.bonds, j)
    
    # Detach from atom graph
    #  Remove at2 from at1.children and set at2.parent to nothing
    # This assumes at1 is parent of at2
    popparent!(at2)

    at1.container === at2.container && return pose

    # Case this is an inter-residue connection, the downstream residue needs to
    # be coupled with the origin and detached from residue graph
    # This assumes ROOT is always on the side of the parent
    # This assumes at1 is parent of at2
    #  Remove at2 from at1.children and set at2.parent to nothing
    #  Add at2 to origin.children and set at2.parent to origin
    _origin = ProtoSyn.root(at1)

    #  Remove at2.container from at1.containter.children and set
    # at2.container.parent to nothing
    hasparent(at2.container) && popparent!(at2.container)
    setparent!(at2, _origin)
    setparent!(at2.container, _origin.container)
    
    # Reindex to set correct ascendents
    reindex(pose.graph)

    return pose
end


"""
    bond(at1::Atom, at2::Atom)
    
Bond both given [`Atom`](@ref) instances (adds `at2` to `at1.bonds` and
vice-versa). Both [`Atom`](@ref) instances need to be in the same
[`Segment`](@ref).

# See also
[`join`](@ref) [`unbond`](@ref)

# Examples
```jldoctest
julia> bond(atom1, atom2)
```
"""
@inline function bond(at1::Atom, at2::Atom)
    @assert at1.container.container === at2.container.container "can only bond atoms within the same segment"
    !in(at2, at1.bonds) && push!(at1.bonds, at2)
    !in(at1, at2.bonds) && push!(at2.bonds, at1)
end


export join
"""
    join(at1::Atom, at2::Atom)

Join [`Atom`](@ref) `at1` with [`Atom`](@ref) `at2`.

    join(r1::Residue, s1::String, r2::Residue, s2::String)
    
Join [`Atom`](@ref) named `s1` from [`Residue`](@ref) `r1` with [`Atom`](@ref)
named `s2` from [`Residue`](@ref) `r2`.

Bond (add eachother to `other.bonds` field) and set parent/children relationship
of both the [`Atom`](@ref) instances and respective `atom.container`
([`Residue`](@ref)). Note that `at2` [`Atom`](@ref) will become parent at `at1`
(and `at2.container` [`Residue`](@ref) will become parent of `at1.container`).

# See also
[`bond`](@ref) [`unbond`](@ref)

# Examples
```jldoctest
julia> join(r1, "C", r2, "N")
```
"""
function join(r1::Residue, s1::String, r2::Residue, s2::String) # ! IMPORTANT
    hasparent(r2) && error("r2 is already connected.")
    at1 = r1[s1]
    at2 = r2[s2]
    bond(at1, at2)          # at1 <-> at2
    setparent!(at2, at1)    # at1 -> at2
    setparent!(r2, r1)      # r1 -> r2
end

function join(at1::Atom, at2::Atom)
    hasparent(at2) && error("at2 is already connected.")
    bond(at1, at2)
    setparent!(at2, at1)
    setparent!(at2.container, at1.container)
end

# COUNTERS
count_segments(t::Topology) = length(t.items)

count_residues(c::AbstractContainer) = mapreduce(x -> count_residues(x), +, c.items; init=0)
count_residues(s::Segment) = length(s.items)
count_residues(r::Residue) = 1

count_atoms(c::AbstractContainer) = mapreduce(x -> count_atoms(x), +, c.items, init=0)
count_atoms(r::Residue) = r.size
count_atoms(a::Atom) = 1


export travel_graph

"""
    travel_graph(start::Atom; [stop::Opt{Atom} = nothing])

Return a `Vector{Atom}` with all atom instances between [`Atom`](@ref) `start`
and `stop`, while following the structure's [Graph](@ref state-types). If no `stop`
[`Atom`](@ref) instance is provided or if it isn't found as a downstream parent
of the `start` [`Atom`](@ref), all instances until no children [`Atom`](@ref)
instances are found are returned (for example, until the end of the current
[Pose](@ref) of [`Segment`](@ref)). Note that the order of the returned
[`Atom`](@ref) instances reflects the organization of the graph followed, and
not the distance/parenthood to the `start` [`Atom`](@ref), and should therefore
be ignored in most cases.

# See also
[`is_contiguous`](@ref) [`hasparent`](@ref) [`setparent!`](@ref)
 
# Examples
```jldoctest
julia> ProtoSyn.travel_graph(pose.graph[1][72][1])
32-element Array{Atom,1}:
 ...
```
"""
function travel_graph(start::Atom, stop::Opt{Atom} = nothing)::Vector{Atom}
    atoms = Vector{Atom}([start])
    stack = Vector{Atom}(copy(start.children))

    while length(stack) > 0
        atom_i = pop!(stack)
        if atom_i != stop
            stack = vcat(stack, copy(atom_i.children))
        end
        push!(atoms, atom_i)
    end

    return atoms
end

export ids

"""
    ids(atoms::Vector{Atom})

Return a vector with the `:id` `Int` field for every [`Atom`](@ref) in the given
`atoms` vector.

# See also
[`travel_graph`](@ref)

# Examples
```jldoctest
ProtoSyn.ids(atoms)
32-element Array{Int64,1}:
 ...
```
"""
function ids(atoms::Vector{Atom})::Vector{Int}
    idxs = Vector{Int}()
    for atom in atoms
        push!(idxs, atom.id)
    end
    return idxs
end


"""
    is_contiguous(pose::Pose, selection::AbstractSelection)

Returns `true` if all the [`Residue`](@ref) instances gathered from the
`selection` applied to the given `pose` are contiguous (have a parenthood
relationship connecting them all). 

# Examples
```jldoctest
julia> ProtoSyn.is_contiguous(pose, rid"1" | rid"3")
false

julia> ProtoSyn.is_contiguous(pose, rid"1:10")
true
```
"""
function is_contiguous(pose::Pose, selection::ProtoSyn.AbstractSelection)
    sele              = promote(selection, Residue)
    selected_residues = sele(pose, gather = true)
    root              = ProtoSyn.root(pose.graph).container # * Residue
    
    # Initialize
    stack = Vector{Residue}()
    for child in root.children
        push!(stack, child)
    end
    
    # Travel graph, search for first selected residue
    selected_stack = Vector{Residue}()
    while length(stack) > 0
        residue = pop!(stack)
        if residue in selected_residues
            push!(selected_stack, residue)
            break
        end
        for child in residue.children
            push!(stack, child)
        end
    end

    length(selected_stack) == 0 && begin
        @warn "Finished travelling the Pose.graph but the selected residues were not found."
        return nothing
    end

    # Travel graph, mark selected atoms found
    marks = 0
    while length(selected_stack) > 0
        marks += 1
        residue = pop!(selected_stack)
        for child in residue.children
            if child in selected_residues
                push!(selected_stack, child)
            end
        end
    end

    # Compare number of found marks with the number of selected residues
    return marks == length(selected_residues)
end