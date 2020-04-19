using LinearAlgebra: dot

export build_tree!

function build_tree!(top::Topology)

    queue = Atom[]
    byatom = eachatom(top)
    
    # reset graph flags
    foreach(r->r.visited=false, eachresidue(top))
    foreach(at->at.visited=false, byatom)

    for atom in byatom
        atom.visited && continue
        atom.visited = true
        push!(queue, atom)
        while !isempty(queue)
            parent = popfirst!(queue)
            for at in parent.bonds
                at.visited && continue
                link!(parent, at)
                at.visited = true
                push!(queue, at)
            end
        end
    end

    root = origin(top)
    for atom in byatom
        
        # if this atom is orphan, then make it a child
        # of the root (origin)
        if !hasparent(atom)
            link!(root, atom)
            continue
        end

        # build residue graph
        child_res = atom.container
        parent_res = atom.parent.container
        if child_res===parent_res || (child_res.visited && parent_res.visited)
            continue
        end
        link!(parent_res, child_res)
        parent_res.visited = true
        child_res.visited = true
    end

    for atom in byatom
        atom.ascendents = ascendents(atom, 4)
    end
    top
end

export ascendents
ascendents(c::AbstractContainer, level::Int) = begin
    level == 1 ? (c.index,) : (c.index, ascendents(c.parent, level-1)...)
end

export sync!
sync!(state::State, top::Topology, force=false) = begin
    if state.c2i && state.i2c
        error("unable to request simultaneous i->c and c->i coordinate conversion")
    elseif state.c2i
        c2i!(state, top)
    elseif state.i2c
        i2c!(state, top, force)
    end
    state
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



i2c!(state::State, top::Topology, force=false) = begin
    # assert top.id==state.id
    
    vjk = MVector{3,Float64}(0.0, 0.0, 0.0)
    vji = MVector{3,Float64}(0.0, 0.0, 0.0)
    n   = MVector{3,Float64}(0.0, 0.0, 0.0)
    
    queue = Atom[]

    root = origin(top)
    xyz = zeros(3, state.size)
    for child in root.children
        # force all child states to be updated
        state[child].changed |= force
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
    # detach from container
    if hascontainer(c)
        delete!(c.container, c)
    end
    
    # detach from graph
    #  1. remove parent
    if hasparent(c)
        unlink!(c.parent, c)
    end
    #  2. remove children
    while !isempty(c.children)
        unlink!(c, pop!(c.children))
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
            deleteat!(atom.bonds, i)
            if (j = findfirst(at->at===atom, other.bonds)) !== nothing
                deleteat!(other.bonds, j)
            end

            # detach from atom graph
            if atom === other.parent
                unlink!(atom, other)
            elseif other === atom.parent
                unlink!(other, atom)
            end
        end

        # remove connection to the root node
        if orig!==nothing && atom.parent===orig
            unlink!(orig, atom)
        end
    end
    r
end


Base.pop!(top::Topology, state::State, res::Residue) = begin
    if state.id != top.id
        error("mismatch between state and topology IDs")
    #elseif isorphan(res)
    #    error("unable to pop orphan residues from topology+state")
    elseif res.container.container !== top
        error("given residue does not belong to the provided topology")
    end

    # detach residue from parents
    detach(res)

    # remove node states from parent state and create
    # a new state for this residue
    st = splice!(state, res[1].index:res[end].index)
    reindex(top)
    
    # new common ID
    res.id = st.id = genid()

    return (res,st)
    
end





#---------------------------
# export append!
Base.append!(segment::Segment, state::State, seq::Vector{String}, db::ResidueDB, rxtb::ReactionToolbelt) = begin
    residues = _insert!(segment, 1, state, 1, seq, db, rxtb)

    root = ProtoSyn.origin(segment)
    link!(root, rxtb.root(residues[1]))

    segment
end

Base.append!(parent::Residue, state::State, seq::Vector{String}, db::ResidueDB, rxtb::ReactionToolbelt) = begin
    segment = parent.container
    
    # residue insertion point
    ripoint = findfirst(r->r===parent, segment.items)
    
    # state insertion point
    sipoint = mapreduce(a->a.index, max, parent.items)
    
    # perform insertion
    residues = _insert!(segment, ripoint+1, state, sipoint+1, seq, db, rxtb)
    
    # link new sequence to parent
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