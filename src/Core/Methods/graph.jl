# This file should contain functions that work on the system Graph, such as 
# functions that deal with parent/child relations, indexes and bonds, among
# others.

export hasparent, haschildren, isparent
export setparent!, popparent!
export popchild!
export root


"""
    root(c::AbstractContainer)

Return the first `Atom` in `AbstractContainer` `c` that has no parent.
"""
root(c::AbstractContainer) = begin
    for atom in eachatom(c)
        !hasparent(atom) && return atom
    end
    nothing
end


"""
    hasparent(c::AbstractDigraph) -> Bool

Test whether the given AbstractDigraph `c` has a parent.
"""
@inline hasparent(c::AbstractDigraph) = c.parent !== nothing


"""
    haschildren(c::AbstractDigraph) -> Bool

Test whether the given AbstractDigraph `c` has children.
"""
@inline haschildren(c::AbstractDigraph) = !isempty(c.children)


"""
    isparent(parent::AbstractDigraph, child::AbstractDigraph) -> Bool

Test whether `parent` is the parent of `child`.
"""
@inline isparent(parent::AbstractDigraph, child::AbstractDigraph) = begin
    parent === child.parent
end

@inline isparent(::Nothing, child::AbstractDigraph) = false


"""
    setparent!(child::T, parent::T) where {T<:AbstractDigraph}

Set `parent` as the parent of `child`, while adding `child` to
`parent.children`.
"""
function setparent!(child::T, parent::T) where {T <: AbstractDigraph}
    hasparent(child) && error("unable to setparent! of non-orphan item")
    push!(parent.children, child)
    child.parent = parent
    child
end


"""
    popparent!(child::AbstractDigraph}

Remove the parent from `child` (sets it to `nothing`) while removing `child`
from `parent.children` (only if `child` is a child of `parent`).
"""
function popparent!(child::AbstractDigraph)
    if hasparent(child)
        parent = child.parent
        i = findfirst(x -> x === child, parent.children)
        if i !== nothing
            deleteat!(parent.children, i)
            child.parent = nothing
        end
    end
    child
end


"""
    popchild!(parent::AbstractDigraph} -> child

Remove child at index 1 from `parent` and return it. Note: The returned element
is orphan. Return `nothing` if the item has no children.
"""
popchild!(parent::AbstractDigraph) = begin
    isempty(parent.children) ? nothing : popparent!(parent.children[1])
end


export hascontainer
@inline hascontainer(c::AbstractContainer) = c.container !== nothing

export origin
@inline origin(t::Topology) = get(t.root, "OO")
@inline origin(c::AbstractContainer) = hascontainer(c) ? origin(c.container) : nothing

export hasgraph
@inline hasgraph(t::Topology) = !isempty(origin(t).children)

export genid
@inline genid() = Int(rand(UInt16))


export reindex

"""
    reindex(t::Topology; set_ascendents = true)
    
Re-indexes the whole topology, setting both the ID and index of elements
inside the topology to the corresponding relative index in the container.items
which they belong to. If `set_ascendents` is set to true (is, by default), each
`Atom` instance `ascendents` field will be updated to reflect the new indices.

# Examples
```jldoctest
julia> reindex(pose.graph)
```
"""
function reindex(t::Topology; set_ascendents = true)
    aid = rid = sid = 0
    for seg in t.items
        seg.id    = (sid += 1)
        seg.index = sid
        for res in seg.items
            res.id    = (rid += 1)
            res.index = rid
            for atm in res.items
                atm.id    = (aid += 1)
                atm.index = aid
            end
        end
    end

    # update ascendents (not possible before because
    # of possible problems with index assignment)
    if set_ascendents
        for atm in eachatom(t)
            atm.ascendents = ascendents(atm, 4)
        end
    end

    t
end


export ascendents

"""
    ascedents(c::AbstractContainer, level::Int)
    
Return a Tuple containing the N (`level`) previous `parent` instances in the
graph (recursivelly).

# Examples
```jldoctest
julia> ascendents(atom, 4)
(1, 2, 3, 4)
```
"""
function ascendents(c::AbstractContainer, level::Int)
    level > 1 ? (c.index, ascendents(c.parent, level-1)...) : (c.index,)
end


"""
    reindex(s::Segment)
    
Re-indexes the given `Segment` instance.

# Examples
```jldoctest
julia> reindex(pose.graph[1])
```
"""
function reindex(s::Segment)
    aid = rid = 0
    for res in s.items
        res.id = (rid += 1)
        res.index = rid
        for atm in res.items
            atm.index = (aid += 1)
        end
    end
    s
end


"""
    unbond(pose::Pose, at1::Atom, at2::Atom)::Pose
    
Return a Pose instance with both given atoms unbonded (removed from eachother
`bonds` list; pops parenthood - if parent - and, if bond is inter-residue, sets
the downstream residue.parent to be the origin of the system, maintaining the 
same relative position as measured from internal coordinates, after sync! is
called).

# Examples
```jldoctest
julia> unbond(pose, atom1, atom2)
```
"""
function unbond(pose::Pose, at1::Atom, at2::Atom)::Pose
    @assert (at2 in at1.bonds) & (at1 in at2.bonds) "Atoms $at1 and $at2 are not bonded and therefore cannot be unbonded."
    isparent(at1, at2) && return _unbond(pose, at1, at2)
    isparent(at2, at1) && return _unbond(pose, at2, at1)
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
    # be coupled with the origin
            
    state = pose.state

    # Detach from residue graph
    # This assumes ROOT is always on the side of the parent
    # This assumes at1 is parent of at2
    #  Remove at2 from at1.children and set at2.parent to nothing
    #  Add at2 to origin.children and set at2.parent to origin
    _origin = ProtoSyn.origin(at1)

    #  Remove at2.container from at1.containter.children and set
    # at2.container.parent to nothing
    popparent!(at2.container)
    setparent!(at2, _origin)
    setparent!(at2.container, _origin.container)
    
    # Reindex to set correct ascendents
    reindex(pose.graph)


    at_s = state[at2]

    # Preserve current placement (based on cartesian coordinates)
    # sync!(pose)
    # println("CALLED on $at2")
    # at_s.b = ProtoSyn.distance(at_s, state[_origin])
    # at_s.θ = ProtoSyn.angle(at_s, state[_origin], state[_origin.parent])
    # at_s.ϕ = ProtoSyn.dihedral(at_s, state[_origin], state[_origin.parent], state[_origin.parent.parent])

    # ProtoSyn.request_i2c(state; all = true)

    return pose
end


"""
    bond(at1::Atom, at2::Atom)
    
Bond both given Atoms (adds at2 to at1.bonds and vice-versa). Both atoms need
to be in the same segment.

# Examples
```jldoctest
julia> bond(atom1, atom2)
```
"""
@inline function bond(at1::Atom, at2::Atom)
    @assert segment(at1) === segment(at2) "can only bond atoms within the same segment"
    !in(at2, at1.bonds) && push!(at1.bonds, at2)
    !in(at1, at2.bonds) && push!(at2.bonds, at1)
end


"""
    join(r1::Residue, s1::String, r2::Residue, s2::String)
    
Join atom named `s1` from residue `r1` with atom named `s2` from residue `r2`.
Bond (add eachother to other.bonds field) and set parent/children of both the
atoms and respective atom.container (residue).

# Examples
```jldoctest
julia> join(r1, "C", r2, "N")
```
"""
function join(r1::Residue, s1::String, r2::Residue, s2::String) # IMPORTANT
    hasparent(r2) && error("r2 is already connected")
    at1 = r1[s1]
    at2 = r2[s2]
    bond(at1, at2)          # at1 <-> at2
    setparent!(at2, at1)    # at1 -> at2
    setparent!(r2, r1)      # r1 -> r2
end


export build_tree!
function build_tree!(seedfinder::Function, top::Topology)

    head = tail = 0
    natoms = count_atoms(top)
    tree = Vector{Atom}(undef, natoms)
    
    #seeds = seedfinder(top)
    root = origin(top)

    # build atom tree
    # for seed in seeds
    for seg in top.items
        seed = seedfinder(seg)
        seed === nothing && continue
        tree[tail+=1] = seed
        seed.visited = true
        while head < tail
            parent = tree[head+=1]
            for atom in parent.bonds
                atom.visited && continue
                setparent!(atom, parent)
                tree[tail+=1] = atom 
                atom.visited = true
            end
        end
        hasparent(seed) && error("invalid seed encountered")
        setparent!(seed, root)
    end
    
    # build residue graph
    if head > 0
        for atom in tree
            atom.ascendents = ascendents(atom, 4)
            r = atom.container
            p = atom.parent.container
            if r===p || p===top.root || (r.visited && p.visited)
                continue
            end
            if r.container !== p.container
                error("parent and child residue must belong to the same segment")
            end
            setparent!(r, p)
            r.visited = p.visited = true
        end
    end
end

# COUNTERS
count_segments(t::Topology) = length(t.items)

count_residues(c::AbstractContainer) = mapreduce(x -> count_residues(x), +, c.items; init=0)
count_residues(s::Segment) = length(s.items)
count_residues(r::Residue) = 1

count_atoms(c::AbstractContainer) = mapreduce(x -> count_atoms(x), +, c.items, init=0)
count_atoms(r::Residue) = r.size
count_atoms(a::Atom) = 1


# POP

"""
    pop!(pose::Pose{Topology}, atom::Atom)::Pose{Atom}

Pop and return the desired `atom` from the given `pose`. In order to do this,
perform the following actions:
 - Unset parents/children
 - Unbond neighbours
 - Remove from graph
 - Remove from state
 - Set new ascendents
 - Update the container itemsbyname

# Examples
```jldoctest
julia> pop!(pose, pose.graph[1][1][2])
```
"""
function pop_atom!(pose::Pose{Topology}, atom::Atom)::Pose{Atom}

    if atom.container.container.container !== pose.graph
        error("Given Atom does not belong to the provided topology.")
    end

    # Save information to return 
    popped_atom = Atom(atom.name, 1, 1, atom.symbol)

    # Save children and parent of this atom (will be removed in next step)
    children    = copy(atom.children)
    grandparent = atom.parent

    # Unset parents/children and unbond neighbours
    for i = length(atom.bonds):-1:1   # Note the reverse loop
        other = atom.bonds[i]
        ProtoSyn.unbond(pose, atom, other)
    end

    # Using saved children and parent, set all child.parent to be this
    # atom.parent
    for child in children
        # child.parent = grandparent
        child.parent = origin(pose.graph)
    end

    # Remove from graph
    deleteat!(atom.container.items, findfirst(atom, atom.container.items))
    atom.container.size -= 1

    # Remove from state
    popped_state = splice!(pose.state, atom.index)

    # Reindex and set ascendents
    reindex(pose.graph)

    # Update container 'itemsbyname'
    pop!(atom.container.itemsbyname, atom.name)

    # Set common ID
    popped_atom.id = popped_state.id = genid()

    return Pose(popped_atom, popped_state)
end

function pop_residue!(pose::Pose{Topology}, residue::Residue)

    if residue.container.container !== pose.graph
        error("Given Residue does not belong to the provided topology.")
    end

    # Remove internal atoms. Notice the inverse loop. Also removes atom-level
    # parenthood and bonds
    for atom in reverse(residue.items)
        pop_atom!(pose, atom)
    end

    # Remove from container.items
    deleteat!(residue.container.items, findfirst(residue, residue.container.items))
    residue.container.size -= 1
end

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

# Detach

Base.detach(s::Segment) = begin
    root = origin(s)
    for at in root.children
        isparent(root, at) && popparent!(at)
    end
    hascontainer(s) && delete!(s.container, s)
end

# --- MIGHT NEED A NEW PLACE

export get_ani_species

"""
TO DO
"""
function get_ani_species(container::AbstractContainer)
    
    periodic_table = Dict("H" => 1, "C" => 6, "N" => 7, "O" => 8, "S" => 16)

    species = Vector{Int64}()
    for atom in eachatom(container)
        push!(species, periodic_table[atom.symbol])
    end

    return species
end

get_ani_species(pose::Pose) = get_ani_species(pose.graph)