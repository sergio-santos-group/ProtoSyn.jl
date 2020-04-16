using LinearAlgebra: dot

export build_tree!

function build_tree!(top::Topology)
    queue = Atom[]
    byatom = eachatom(top)
    
    # reset graph flags
    foreach(r->r.node.visited=false, eachresidue(top))
    foreach(at->at.node.visited=false, byatom)

    for atom in byatom
        atom.node.visited && continue
        atom.node.visited = true
        push!(queue, atom)
        while !isempty(queue)
            parent = popfirst!(queue)
            for at in parent.bonds
                at.node.visited && continue
                push!(parent.node, at.node)
                at.node.visited = true
                push!(queue, at)
            end
        end
    end

    root = origin(top).node
    for atom in byatom
        node = atom.node

        # if this atom is orphan, then make it a child
        # of the root (origin)
        if !hasparent(node)
            push!(root, node)
            continue
        end

        # build residue graph
        child_res = node.item.container
        parent_res = node.parent.item.container
        if child_res===parent_res || (child_res.node.visited && parent_res.node.visited)
            continue
        end
        push!(parent_res.node, child_res.node)
        parent_res.node.visited = true
        child_res.node.visited = true
    end

    for atom in byatom
        atom.node.ascendents = ascendents(atom.node, 4)
        # node = atom.node
        # node.ascendents = (
        #     node.item.index,
        #     node.parent.item.index,
        #     node.parent.parent.item.index,
        #     node.parent.parent.parent.item.index
        # )
    end
    top
end

export ascendents
ascendents(n::GraphNode, level::Int) = begin
    if level == 1
        return (n.item.index,)
    end
    (n.item.index, ascendents(n.parent, level-1)...)
end

export sync!
sync!(state::State, top::Topology, force=false) = begin
    if state.c2i && state.i2c
        throw(Exception)
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
        (i,j,k,l) = atom.node.ascendents
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
    
    queue = AtomGraphNode[]

    xyz = zeros(3, state.size)
    root = origin(top).node
    # append!(queue, root.children)
    for child in root.children
        state[child.item.index].changed |= force
        push!(queue, child)
    end
    
    while !isempty(queue)
        node = popfirst!(queue)
        (i,j,k) = node.ascendents
        # i = node.item.index
        istate = state[i]
        println(node)
        for child in node.children
            state[child.item.index].changed |= istate.changed
            push!(queue, child)
        end
        !(istate.changed) && continue
        istate.changed = false
        #println("updating node $i")

        # j = node.parent.item.index
        # k = node.parent.parent.item.index
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



#----------------------
# PEPTIDES PEPTIDES PEPTIDES PEPTIDES PEPTIDES PEPTIDES
#rootprovider(r::Residue) = get(r, "N")

#const ResidueDB = Dict{String, Tuple{Residue, State}}

# function load(dir::AbstractString)
    
#     lib = ResidueDB()

#     files = filter(f->endswith(f, ".pdb"), readdir(dir))
#     foreach(files) do fname
#         top,state = read(joinpath(dir,fname), PDB)
#         foreach(eachresidue(top)) do residue
#             lib[residue.name] = pop!(top, state, residue)
#         end
#     end
#     lib
# end

# # Base.pop!(seg::Segment, res::Residue) = begin
# #     i = findfirst(segres.segment
# # end
# Base.findfirst(ar::Vector{T}, item::T) where T = begin
# Base.findfirst(item::T, ar::Vector{T}) where T = begin
#     return findfirst(it->it===item, ar)
# end

@inline Base.in(item::T, container::AbstractContainer{T}) where T = begin
    item.container===container
end

Base.findfirst(item::T, container::AbstractContainer{T}) where T = begin
    in(item,container) ? findfirst(x->x===item, container.items) : nothing
end
Base.findfirst(item::GraphNode{T}, container::GraphNode{T}) where T = 
    findfirst(x->x===item, container.children)

# @inline Base.in(a::Atom, r::Residue) = a.container===r
# @inline Base.in(r::Residue, s::Segment) = r.container===s

# Base.findfirst(a::Atom, r::Residue) = begin
#     in(a,r) ? findfirst(atm->atm===a, r.items) : nothing
# end

# Base.findfirst(r::Residue, s::Segment) = begin
#     in(r,s) ? findfirst(res->res===r, s.residues) : nothing
# end

Base.delete!(container::AbstractContainer{T}, item::T) where T = begin
    if (i = findfirst(item, container)) !== nothing
        deleteat!(container.items, i)
        item.container = nothing
        container.size -= 1
    end
    container
end

Base.delete!(parent::GraphNode{T}, child::GraphNode{T}) where T = begin
    if (i = findfirst(child, parent))!== nothing
        deleteat!(parent.children, i)
        child.parent = nothing
    end
    parent
end

Base.detach(node::GraphNode{T}) where T = begin
    hasparent(node) && delete!(node.parent, node)
    while !isempty(node.children)
        pop!(node.children).parent=nothing
    end
    node
end

# intop(c::AbstractContainer) = hasparent(c) && hasparent(c.container)
# intop(t::Topology) = true

Base.detach(r::Residue) = begin
    orig = origin(r)

    # detach from segment
    if hasparent(r)
        delete!(r.container, r)
    end
    
    # detach from residue graph
    detach(r.node)
    
    # remove all inter-residue atom bonds
    for atom in eachatom(r)
        node = atom.node
        for i=length(atom.bonds):-1:1
            other = atom.bonds[i]
            # do nothing if the other atom is in this residue
            in(other,r) && continue

            # remove inter-residue bonds
            deleteat!(atom.bonds, i)
            if (j = findfirst(at->at===atom, other.bonds)) !== nothing
                deleteat!(other.bonds, j)
            end
            # detach from atom graph
            if node === other.node.parent
                delete!(node, other.node)
            elseif other.node === atom.node.parent
                delete!(other.node, node)
            end
        end
        if orig!==nothing && atom.node.parent===orig.node
            delete!(orig.node, atom.node)
        end
    end
    r
end

Base.parent(c::AbstractContainer{T}) where T = c.node.parent

Base.pop!(top::Topology, state::State, res::Residue) = begin
    if state.id != top.id
        error("mismatch between state and topology IDs")
    elseif isorphan(res)
        error("unable to pop orphan residues from topology+state")
    elseif res.container.container !== top
        error("given residue does not belong to the provided topology")
    end

    # detach residue from parents
    detach(res)

    # remove node states from parent state and create
    # a new state for this residue
    atoms = res.items
    st = splice!(state, atoms[1].index:atoms[end].index)
    
    # new common ID
    res.id = st.id = genid()

    # renumber
    reindex(top)

    return (res,st)
    
end


# # Base.Colon(a::Atom, b::Atom) = UnitRange(a.index,b.index)
# # import Base.-

# # :(a::AbstractContainer{T}, b::AbstractContainer{T}) where T = a.index:b.index
# Base.:-(a::AbstractContainer{T}, b::AbstractContainer{T}) where T = a.index-b.index
# # Base.:-(a::Atom, b::Atom) = a.index-b.index
# Base.Colon(a::Atom, b::Atom) = UnitRange{Int}(a.index,b.index)

Base.copy(rs::Tuple{Residue, State}) = (copy(rs[1]), copy(rs[2]))