using LinearAlgebra: dot

export build_tree!

export bond, unbond

# residue(a::Atom) = a.container
# segment(r::Residue) = r.container
# topology(s::Segment) = s.container

segment(at::Atom) = hascontainer(at) ? at.container.container : nothing

Base.findfirst(x::T, c::Vector{T}) where T = findfirst(i->i===x, c)
Base.findfirst(x::T, c::AbstractContainer{T}) where T = findfirst(x, c.items)

@inline function unbond(at1::Atom, at2::Atom)
    i = findfirst(at1, at2.bonds)
    i !== nothing && deleteat!(at2.bonds, i)
    
    j = findfirst(at2, at1.bonds)
    j !== nothing && deleteat!(at1.bonds, j)
end

@inline function bond(at1::Atom, at2::Atom)
    @assert segment(at1)===segment(at2) "can only bond atoms within the same segment"
    !in(at2, at1.bonds) && push!(at1.bonds, at2)
    !in(at1, at2.bonds) && push!(at2.bonds, at1)
end



# function build_tree(molecule::ProtoSyn.Molecule)
#     ts = ProtoSyn.TraverseState(size(molecule))
#     head = tail = 0

#     root = TreeNode(get(molecule.residues[1], "N"))
#     tree = Vector{TreeNode}(undef, size(molecule))
    
#     ts.visited[root.atom.id+root.atom.residue.offset] = true
#     tree[tail+=1] = root

#     #root.id+root.residue.offset
    
#     atoms = collect(ProtoSyn.iterbyatom(molecule))

#     while head < tail
#         # id1 = ts.indices[head+=1]
#         parent = tree[head+=1]
#         id1 = parent.atom.id + parent.atom.residue.offset
#         for id2 in molecule.bonds[id1]
#             if !ts.visited[id2]
#                 node = TreeNode(atoms[id2], root)
#                 push!(parent.children, node)
#                 ts.visited[id2] = true
#                 # ts.indices[tail+=1] = id2
#                 tree[tail+=1] = node
#             end
#         end
#     end
#     return tree
# end

#build_tree!(top::Topology) = build_tree!(t->[t[1,1,"N"]], top)

function build_tree!(seedfinder::Function, top::Topology)

    head = tail = 0
    natoms = count_atoms(top)
    tree = Vector{Atom}(undef, natoms)
    
    seeds = seedfinder(top)
    root = origin(top)

    # build atom tree
    for seed in seeds    
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
    tree
end

count_atoms(c::AbstractContainer) = mapreduce(x->count_atoms(x), +, c.items, init=0)
count_atoms(r::Residue) = r.size
count_atoms(a::Atom) = 1


# function build_tree!(top::Topology)

#     queue = Atom[]
#     byatom = eachatom(top)
    
#     # reset graph flags
#     foreach(r->r.visited=false, eachresidue(top))
#     foreach(at->at.visited=false, byatom)

#     for atom in byatom
#         atom.visited && continue
#         atom.visited = true
#         push!(queue, atom)
#         while !isempty(queue)
#             parent = popfirst!(queue)
#             for at in parent.bonds
#                 at.visited && continue
#                 setparent!(at, parent)
#                 at.visited = true
#                 push!(queue, at)
#             end
#         end
#     end

#     root = origin(top)
#     for atom in byatom
        
#         # if this atom is orphan, then make it a child
#         # of the root (origin)
#         if !hasparent(atom)
#             setparent!(atom, root)
#             continue
#         end

#         # build residue graph
#         child_res = atom.container
#         parent_res = atom.parent.container
#         if child_res===parent_res || (child_res.visited && parent_res.visited)
#             continue
#         end
#         setparent!(child_res, parent_res)
#         parent_res.visited = true
#         child_res.visited = true
#     end

#     for atom in byatom
#         atom.ascendents = ascendents(atom, 4)
#     end
#     top
# end


export ascendents
ascendents(c::AbstractContainer, level::Int) = begin
    level > 1 ? (c.index, ascendents(c.parent, level-1)...) : (c.index,)
    # level == 1 ? (c.index,) : (c.index, ascendents(c.parent, level-1)...)
end

export sync!
# sync!(state::State, top::Topology) = begin
sync!(pose::Pose{Topology}) = begin
    state = pose.state
    top = pose.graph
    if state.c2i && state.i2c
        error("unable to request simultaneous i->c and c->i coordinate conversion")
    elseif state.c2i
        c2i!(state, top)
    elseif state.i2c
        i2c!(state, top)
    end
    pose
end





c2i!(state::State{T}, top::Topology) where T = begin
    vij = MVector{3,T}(T(0), T(0), T(0))
    vjk = MVector{3,T}(T(0), T(0), T(0))
    vkl = MVector{3,T}(T(0), T(0), T(0))
    n   = MVector{3,T}(T(0), T(0), T(0))
    m   = MVector{3,T}(T(0), T(0), T(0))
    o   = MVector{3,T}(T(0), T(0), T(0))
    for atom in eachatom(top)
        (i,j,k,l) = atom.ascendents
        istate = state[i]
        
        # bond
        @. vij = state[j].t - istate.t
        dij = sqrt(dot(vij,vij))
        istate.b = dij

        # angle
        @. vjk = state[k].t - state[j].t
        djk = sqrt(dot(vjk,vjk))
        istate.θ = pi - acos(dot(vij,vjk) / (dij*djk))

        # dihedral
        @. vkl = state[l].t - state[k].t
        @cross u n[u] vij[u] vjk[u]
        @cross u m[u] vjk[u] vkl[u]
        @cross u o[u] n[u] m[u]
        x = dot(o,vjk)/sqrt(dot(vjk,vjk))
        y = dot(n,m)
        istate.ϕ = atan(x,y)
    end
    state.c2i = false
    state
end



i2c!(state::State, top::Topology) = begin
    # assert top.id==state.id
    
    vjk = MVector{3,Float64}(0.0, 0.0, 0.0)
    vji = MVector{3,Float64}(0.0, 0.0, 0.0)
    n   = MVector{3,Float64}(0.0, 0.0, 0.0)
    
    queue = Atom[]

    root = origin(top)
    root_changed = state[root].changed
    xyz = zeros(3, state.size)
    for child in root.children
        # force all child states to be updated
        state[child].changed |= root_changed
        push!(queue, child)
    end
    
    while !isempty(queue)
        atom = popfirst!(queue)
        (i,j,k) = atom.ascendents
        
        istate = state[i]
        for child in atom.children
            state[child].changed |= istate.changed
            push!(queue, child)
        end
        !(istate.changed) && continue
        istate.changed = false
        
        jstate = state[j]        
        kstate = state[k]        
        Ri = istate.r
        Rj = jstate.r
        
        # local coord system
        b = istate.b
        sθ,cθ = sincos(istate.θ)  # angle
        sϕ,cϕ = sincos(istate.ϕ + jstate.Δϕ)  # dihedral
        x_1 = -b*cθ
        x_2 =  b*cϕ*sθ
        x_3 =  b*sϕ*sθ
        
        # rotate to parent coord system
        @nexprs 3 u -> vji[u] = Rj[u,1]*x_1 + Rj[u,2]*x_2 + Rj[u,3]*x_3
        
        # UPDATE ROTATION MATRIX
        # @nexprs 3 u -> vjk_u = kstate.t[u] - jstate.t[u]
        @. vjk = kstate.t - jstate.t
        
        # column 1 (x)
        @nexprs 3 u -> Ri[u,1] = vji[u]/b
            
        # column 3 (z)
        @cross u n[u] vji[u] vjk[u]
        dn = sqrt(dot(n,n))
        @nexprs 3 u -> Ri[u,3] = n[u]/dn
    
        # column 2 (y)
        @cross u Ri[u,2] Ri[u,3] Ri[u,1]
        
        # move to new position
        # @nexprs 3 u -> istate.t[u] = vji_u + jstate.t[u]
        # xyz[:, i] .= istate.t
        @. istate.t = vji + jstate.t
        @. xyz[:,i] = istate.t
    end
    xyz
end


function _detach(c::AbstractContainer)
    # 1. detach from container
    hascontainer(c) && delete!(c.container, c)
    
    # 2. detach from graph
    #  2.1. remove parent
    hasparent(c) && popparent!(c)
    
    #  2.2 remove children
    while haschildren(c)
        popchild!(c)
    end
    c
end


Base.detach(r::Residue) = begin
    # identify the origin of this residue before detachment,
    # otherwise it will no longer be possible!
    orig = origin(r)

    _detach(r)

    # remove all inter-residue atom bonds
    for atom in eachatom(r)
        for i=length(atom.bonds):-1:1   # note the reverse loop
            other = atom.bonds[i]

            # do nothing if the other atom is in this residue
            in(other,r) && continue

            # remove inter-residue bonds
            unbond(atom, other)

            # detach from atom graph
            # if atom === other.parent
            if isparent(atom, other)
                popparent!(other)
            # elseif other === atom.parent
            elseif isparent(other, atom)
                popparent!(atom)
            end
        end

        # remove connection to the root node
        if orig!==nothing && isparent(orig, atom)
            popparent!(atom)
        end
    end
    r
end


# Base.pop!(top::Topology, state::State, res::Residue) = begin
Base.pop!(pose::Pose{Topology}, res::Residue) = begin
    #if state.id != top.id
    #    error("mismatch between state and topology IDs")
    #elseif isorphan(res)
    #    error("unable to pop orphan residues from topology+state")
    if res.container.container !== pose.graph
        error("given residue does not belong to the provided topology")
    end

    # detach residue from parents
    detach(res)

    # remove node states from parent state and create
    # a new state for this residue
    st = splice!(pose.state, res[1].index:res[end].index)
    reindex(pose.graph)
    
    # new common ID
    res.id = st.id = genid()

    #return (res,st)
    Pose(res,st)
end


function fragment(pose::Pose{Topology})
    
    length(pose.graph) != 1 && error("only topologies with a single segment can be turned into fragments")
    
    topology = pose.graph
    segment = topology[1]
    (imin,imax) = extrema(map(at->at.index, eachatom(segment)))
    state = splice!(pose.state, imin:imax)
    detach(segment)
    segment.id = state.id = genid()

    Pose(segment, state)
end
Base.detach(s::Segment) = hascontainer(s) && delete!(s.container, s)


# function fragment(pose::Pose{Segment})
    
#     length(pose.graph) != 1 && error("only segments with a single residue can be turned into fragments")
    
#     segment = pose.graph
#     residue = segment[1]
#     state = splice!(pose.state, residue)
#     detach(residue)
#     residue.id = state.id = genid()

#     Pose(residue, state)
# end




# const Fragment = Vector{Residue}
# isfragment(p::Pose) = !hascontainer(p.graph) && !isempty(p.graph)
isfragment(p::Pose) = !(hascontainer(p.graph) || isempty(p.graph))


#---------------------------
"""
append a fragment as a new segment
"""
Base.append!(pose::Pose{Topology}, frag::Pose{Segment}, rxtb::ReactionToolbelt) = begin
    !isfragment(frag) && error("invalid fragment")
    push!(pose.graph, frag.graph)
    append!(pose.state, frag.state)
    
    setparent!(
        rxtb.root(frag.graph),
        origin(pose.graph)
    )
    reindex(pose.graph)
    pose
end

"""
append a fragment onto an existing segment
"""
# Base.append!(pose::Pose{Segment}, frag::Pose{Segment}, rxtb::ReactionToolbelt) = begin

#     !isfragment(frag) && error("invalid fragment")
    
#     # state insertion point. If this segment is empty,
#     # then go to all previous segments to identify the first
#     # non-empty one.
#     segment = pose.graph
#     topology = segment.container
#     i = findfirst(segment, topology)
#     while i > 1 && isempty(topology[i])
#         i -= 1
#     end
#     sipoint = mapreduce(a->a.index, max, eachatom(topology[i]); init=0)
#     insert!(pose.state, sipoint+1, frag.state)

#     # simply append residues to this segment
#     residues = frag.graph.items
#     push!(segment, residues...)
    
#     root = origin(pose.graph)
#     setparent!(rxtb.root(frag.graph), root)
#     reindex(pose.graph)
#     pose
# end







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

Base.append!(parent::Residue, state::State, seq::Vector{String}, db::ResidueDB, rxtb::ReactionToolbelt) = begin
    segment = parent.container
    
    # residue insertion point
    ripoint = findfirst(parent, segment)
    # ripoint = findfirst(r->r===parent, segment.items)
    
    # state insertion point
    sipoint = mapreduce(a->a.index, max, parent.items)
    
    # perform insertion
    residues = _insert!(segment, ripoint+1, state, sipoint+1, seq, db, rxtb)
    
    # join new sequence to parent
    rxtb.join(parent, residues[1])
    segment
end


_insert!(segment::Segment, ripoint::Int, state::State, sipoint::Int, seq::Vector{String}, db::ResidueDB, rxtb::ReactionToolbelt) = begin
    prev = nothing
    residues = Residue[]
    for s in seq
        res, st = from(db, s)
        
        # insert sub-state into state
        insert!(state, sipoint, st)
        sipoint += length(res)
        
        # insert residue into segment
        insert!(segment, ripoint, res)
        res.id = ripoint
        ripoint += 1

        prev !== nothing && rxtb.join(prev, res)
        push!(residues, res)
        prev = res
    end
    residues
end

Base.insert!(segment::Segment, state::State, i::Int, j::Int, seq::Vector{String}, db::ResidueDB, rxtb::ReactionToolbelt) = begin
    if i < 0
        throw(BoundsError(segment, i))
    elseif i > length(segment)
        throw(BoundsError(segment, j))
    end

    r1 = segment[i]
    r2 = segment[j]
    
    sipoint = mapreduce(a->a.index, max, r1.items; init=0)
    ripoint = findfirst(r->r===r1, segment.items)
    residues = _insert!(segment, ripoint+1, state, sipoint+1, seq, db, rxtb)
    rxtb.join(r1, residues[1])
    rxtb.join(residues[end], r2)
    
    segment
end




export setdihedral!

@inline setdihedral!(s::State{T}, r::Residue, atname::AbstractString, value::T) where T = begin
    at = r[atname]
    at !== nothing && setdihedral!(s, at, value)
    s
end

@inline setdihedral!(s::State{T}, at::Atom, val::T) where T = (s[at].Δϕ = val; s)


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
