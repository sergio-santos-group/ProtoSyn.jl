


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


function join(r1::Residue, s1::String, r2::Residue, s2::String) # IMPORTANT
    hasparent(r2) && error("r2 is already connected")
    at1 = r1[s1]
    at2 = r2[s2]
    bond(at1, at2)          # at1 <-> at2
    setparent!(at2, at1)    # at1 -> at2
    setparent!(r2, r1)      # r1 -> r2
end