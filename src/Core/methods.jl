


function Base.pop!(pose::Pose{Topology}, atom::Atom)::Pose{Atom}
    # Should:
    # - Unset parents/children
    # - Unbond neighbours
    # - Remove from graph
    # - Remove from state
    # - Set new ascendents
    # - Update the container itemsbyname

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
        unbond(pose, atom, other)
    end

    # Using saved children and parent, set all child.parent to be this
    # atom.parent
    for child in children
        child.parent = grandparent
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


function Base.pop!(pose::Pose{Topology}, residue::Residue)

    for atom in residue.items
        pop!(pose, atom)
    end
end


# Base.pop!(pose::Pose{Topology}, r::Residue) = begin

#     println("Deleting residue $r")

#     if r.container.container !== pose.graph
#         error("given residue does not belong to the provided topology")
#     end

#     # detach residue from parents (bonds, parents/children)
#     # hascontainer(r) && delete!(r.container, r)
#     for atom in eachatom(r)
#         for i = length(atom.bonds):-1:1   # Note the reverse loop
#             other = atom.bonds[i]

#             in(other, r) && continue
#             Builder.unbond(pose, atom, other)
#         end
#     end

#     # remove node states from parent state and create
#     # a new state for this residue
#     deleteat!(r.container.items, findfirst(r, r.container.items))
#     r.container.size -= 1
#     for child in origin(pose.graph).children
#         child in r.items && popparent!(child) 
#     end

#     st = splice!(pose.state, r[1].index:r[end].index)
#     reindex(pose.graph) # Also sets ascendents
    
#     # new common ID
#     r.id = st.id = genid()

#     #return (res, st)
#     Pose(r, st)
# end

Base.pop!(pose::Pose{Topology}, seg::Segment) = begin

    if seg.container !== pose.graph
        error("given residue does not belong to the provided topology")
    end

    # remove node states from parent state and create
    # a new state for this segment
    st = splice!(pose.state, seg[1][1].index:seg[end][end].index)
    
    pop!(pose.graph.items, seg)
    # deleteat!(pose.graph.items, findall(x -> x == seg, pose.graph.items))
    # pose.graph.size -= 1

    reindex(pose.graph)
    
    # new common ID
    seg.id = st.id = genid()
    Pose(seg, st)
end


Base.pop!(top::Topology, seg::Segment) = begin
    deleteat!(top.items, findall(x -> x == seg, top.items))
    top.size -= 1
end


function fragment(pose::Pose{Topology})
    
    length(pose.graph) != 1 && error("only topologies with a single segment can be turned into fragments")
    
    topology = pose.graph
    segment = topology[1]
    #(imin,imax) = extrema(map(at->at.index, eachatom(segment)))
    #state = splice!(pose.state, imin:imax)
    state = splice!(pose.state, 1:count_atoms(segment))
    detach(segment)
    segment.id = state.id = genid()
    segment.name = topology.name

    Pose(segment, state)
end

Base.detach(s::Segment) = begin
    root = origin(s)
    for at in root.children
        isparent(root, at) && popparent!(at)
    end
    hascontainer(s) && delete!(s.container, s)
end

root(c::AbstractContainer) = begin
    for atom in eachatom(c)
        !hasparent(atom) && return atom
    end
    nothing
end

isfragment(p::Pose) = !(hascontainer(p.graph) || isempty(p.graph))

export append
function append end


#---------------------------
"""
    append!(pose::Pose{Topology}, frag::Fragment)

Append a fragment as a new segment.
"""
Base.append!(pose::Pose{Topology}, frag::Fragment) = begin
    # NOTE: THIS IS THE APPEND CALLED BY BUILD.
    # Note: A Fragment is a Pose with just 1 Segment.

    !isfragment(frag) && error("invalid fragment")
    
    # Merge the fragment graph (Segment) to the pose graph (Topology).
    push!(pose.graph, frag.graph)

    # Merge the fragment state to the pose state.
    Base.append!(pose.state, frag.state)
    
    # Make sure the fragment graph has the same origin of the new pose.
    root_residue = root(frag.graph).container
    setparent!(
        root(frag.graph),
        origin(pose.graph)
    )

    setparent!(
        root_residue,
        origin(pose.graph).container
    )

    # Re-index the pose to account for the new segment/residue/atoms
    reindex(pose.graph)
    pose
end

# append(pose::Pose{Topology}, frag::Fragment, rxtb::ReactionToolbelt) = begin
#     !isfragment(frag) && error("invalid fragment")
    
#     push!(pose.graph, frag.graph)
#     append!(pose.state, frag.state)
    
#     setparent!(
#         rxtb.root(frag.graph),
#         origin(pose.graph)
#     )
#     reindex(pose.graph)
#     pose
# end


"""
    append!(pose::Pose{Topology}, frag::Fragment, segment::Segment, rxtb::ReactionToolbelt) = begin

Append a fragment onto an existing segment. This fragment will be made a child
to the origin of the provided topology.
"""
append!(pose::Pose{Topology}, frag::Fragment, segment::Segment, rxtb::ReactionToolbelt) = begin
    
    if !in(segment, pose.graph)
        error("the given segment must belong to the given topology")
    elseif !isfragment(frag)
        error("invalid fragment")
    end

    sip = sipoint(segment)
    _insert(segment, pose.state, frag, sip)
    
    setparent!(
        rxtb.root(frag.graph),
        origin(pose.graph)
    )

    reindex(pose.graph)
    pose
end


"""
    append!(pose::Pose{Topology}, frag::Fragment, residue::Residue, rxtb::ReactionToolbelt)

Append a `fragment` to the containing segment of the provided `residue` and
make it a child of that same `residue`.
"""
append!(pose::Pose{Topology}, frag::Fragment, residue::Residue, rxtb::ReactionToolbelt) = begin
        
    if !in(residue, pose.graph)
        error("the given residue must belong to the given topology")
    elseif !isfragment(frag)
        error("invalid fragment")
    end

    # state insertion point
    segment = residue.container
    sip = sipoint(segment)

    _insert(segment, pose.state, frag, sip)
    rxtb.join(residue, frag.graph[1])
    
    reindex(pose.graph)
    pose
end

export insert
insert(pose::Pose{Topology}, frag::Fragment, r1::Residue, r2::Residue, rxtb::ReactionToolbelt) = begin
    if !in(r1, pose.graph)
        error("r1 must belong to the given topology")
    elseif !in(r2, pose.graph)
        error("r2 must belong to the given topology")
    elseif r1.container !== r2.container
        error("r1 and r2 must be in the same segment")
    elseif r1 !== r2.parent
        error("r2 must be a child of r1")
    elseif !isfragment(frag)
        error("invalid fragment")
    end
    
    segment = r1.container
    sip = sipoint(segment)
    rxtb.split(r1, r2)
    _insert(segment, pose.state, frag, sip)
    rxtb.join(r1, frag.graph[1])
    rxtb.join(r2, frag.graph[end])
    
    reindex(pose.graph)
    pose
end

"""
    sipoint(seg::Segment)

Determine the insertion point for the given `segment`.
"""
function sipoint(seg::Segment)
    top = seg.container
    i = findfirst(seg, top)
    while i > 1 && isempty(top[i])
        i -= 1
    end
    lstidx = isempty(top[1]) ? 0 : mapreduce(a->a.index, eachatom(top[1,end]))
    # lstidx = mapreduce(a->a.index, max, eachatom(top[i]); init=0)
    lstidx+1
end


# append(pose::Pose{Topology}, frag::Fragment, residue::Residue, rxtb::ReactionToolbelt) = begin
        
#     if !in(residue, pose.graph)
#         error("the given anchor must belong to the given topology")
#     elseif !isfragment(frag)
#         error("invalid fragment")
#     end
    
#     # state insertion point
#     segment = residue.container
#     sip = sipoint(segment)

#     _insert(segment, pose.state, frag, sip)
#     rxtb.join(residue, frag.graph[1])
    
#     reindex(pose.graph)
#     pose
# end


function _insert(segment::Segment, state::State, frag::Fragment, sip::Int)
    push!(segment, frag.graph.items...)
    insert!(state, sip, frag.state)
end




# export append!
# Base.append!(segment::Segment, state::State, seq::Vector{String}, db::ResidueDB, rxtb::ReactionToolbelt) = begin
#     residues = _insert!(segment, 1, state, 1, seq, db, rxtb)
# Base.append!(segment::Segment, state::State, frag::Fragment, fstate::State, rxtb::ReactionToolbelt) = begin
#     !isfragment(frag) && error("invalid fragment")
    
#     # state insertion point. If this segment is empty,
#     # then go to all previous segments to identify the first
#     # non-empty one.
#     topology = segment.container
#     i = findfirst(segment, topology)
#     while i > 1 && isempty(topology[i])
#         i -= 1
#     end
#     sipoint = mapreduce(a->a.index, max, eachatom(topology[i]); init=0)
#     insert!(state, sipoint+1, fstate)

#     # simply append residues to this segment
#     foreach(r->push!(segment, r), frag)
    
#     root = origin(segment)
#     setparent!(rxtb.root(frag[1]), root)
    
#     segment
# end

# Base.append!(parent::Residue, state::State, seq::Vector{String}, db::ResidueDB, rxtb::ReactionToolbelt) = begin
#     segment = parent.container
    
#     # residue insertion point
#     ripoint = findfirst(parent, segment)
#     # ripoint = findfirst(r->r===parent, segment.items)
    
#     # state insertion point
#     sipoint = mapreduce(a->a.index, max, parent.items)
    
#     # perform insertion
#     residues = _insert!(segment, ripoint+1, state, sipoint+1, seq, db, rxtb)
    
#     # join new sequence to parent
#     rxtb.join(parent, residues[1])
#     segment
# end


# _insert!(segment::Segment, ripoint::Int, state::State, sipoint::Int, seq::Vector{String}, db::ResidueDB, rxtb::ReactionToolbelt) = begin
#     prev = nothing
#     residues = Residue[]
#     for s in seq
#         res, st = from(db, s)
        
#         # insert sub-state into state
#         insert!(state, sipoint, st)
#         sipoint += length(res)
        
#         # insert residue into segment
#         insert!(segment, ripoint, res)
#         res.id = ripoint
#         ripoint += 1

#         prev !== nothing && rxtb.join(prev, res)
#         push!(residues, res)
#         prev = res
#     end
#     residues
# end

# Base.insert!(segment::Segment, state::State, i::Int, j::Int, seq::Vector{String}, db::ResidueDB, rxtb::ReactionToolbelt) = begin
#     if i < 0
#         throw(BoundsError(segment, i))
#     elseif i > length(segment)
#         throw(BoundsError(segment, j))
#     end

#     r1 = segment[i]
#     r2 = segment[j]
    
#     sipoint = mapreduce(a->a.index, max, r1.items; init=0)
#     ripoint = findfirst(r->r===r1, segment.items)
#     residues = _insert!(segment, ripoint+1, state, sipoint+1, seq, db, rxtb)
#     rxtb.join(r1, residues[1])
#     rxtb.join(residues[end], r2)
    
#     segment
# end


export setdihedral!
@inline setdihedral!(s::State, at::Atom, val::T) where {T <: AbstractFloat} = begin
    # Sets dihedral so it becomes the value 'val' 
    # println("Setting dihedral at atom $at to $(rad2deg(val)) (Ascendents: $(at.ascendents))")
    s[at].Δϕ = val - s[at].ϕ
    ProtoSyn.request_i2c(s)
    s
end


export getdihedral!
@inline getdihedral!(s::State, at::Atom) = begin
    return s[at].ϕ + s[at].Δϕ
end


export setoffset!

setoffset!(state::State{T}, at::Atom, default::Number) where T = begin
    # rotates all sibling dihedrals of "at" so that the
    # dihedral angle identified by "at" is equal to "default" 
    if hasparent(at)
        ϕ = state[at].ϕ - T(default)
        for child in at.parent.children
            state[child].ϕ -= ϕ
        end
    end
    state.i2c = true
    state
end


# translate!(state::State, s::Segment, dx::Number) = begin
#     for at in eachatom(s)
#         state[at].t += dx
#     end
#     state.c2i = true
#     state
# end

# rotate!(state::State, s::Segment, r::Matrix) = begin
#     for at in eachatom(s)
#         t = state[at].t
#         t .= r*t
#     end
#     state.c2i = true
#     state
# end

function join(r1::Residue, s1::String, r2::Residue, s2::String) # IMPORTANT
    # println("R1: $r1")
    # println("R2: $r2")
    hasparent(r2) && error("r2 is already connected")
    at1 = r1[s1]
    at2 = r2[s2]
    bond(at1, at2)          # at1 <-> at2
    setparent!(at2, at1)    # at1 -> at2
    setparent!(r2, r1)      # r1 -> r2
end