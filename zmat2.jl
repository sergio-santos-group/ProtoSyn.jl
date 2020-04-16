module Zmat

using StaticArrays
using LinearAlgebra: dot


include("src/Core/macros.jl")

export Atom, Residue, Segment, Topology
export AtomGraphNode, ResidueGraphNode

export Opt
const Opt = Union{Nothing, T} where T

abstract type AbstractIdentifiable end
abstract type AbstractAtom <: AbstractIdentifiable end
abstract type AbstractResidue <: AbstractIdentifiable end
abstract type AbstractSegment <: AbstractIdentifiable end
abstract type AbstractTopology <: AbstractIdentifiable end

mutable struct GraphNode{T<:AbstractIdentifiable}
    item::T
    parent::Opt{GraphNode{T}}
    children::Vector{GraphNode{T}}
    visited::Bool
    ascendents::Opt{NTuple{4,Int}}
end
GraphNode{T}(item::T) where T = GraphNode(item, nothing, Vector{GraphNode{T}}(), false, nothing)

mutable struct Atom <: AbstractAtom
    name::String
    id::Int
    index::Int
    symbol::String
    bonds::Vector{Atom}
    residue::AbstractResidue
    node::GraphNode{Atom}
    Atom(name::String, id::Int, index::Int, symbol::String) = begin
        at = new(name, id, index, symbol, Atom[])
        # at.node = GraphNode{Atom}(at, nothing, [], false)
        at.node = GraphNode{Atom}(at)
        at
    end
end

mutable struct Residue <: AbstractResidue
    name::String
    id::Int
    atoms::Vector{Atom}
    atomsbyname::Dict{String, Atom}
    
    segment::AbstractSegment
    node::GraphNode{Residue}
    Residue(name::String, id::Int) = begin
        r = new(name, id, Atom[], Dict{String,Atom}())
        # r.node = GraphNode{Residue}(r, nothing, [], false)
        r.node = GraphNode{Residue}(r)
        r
    end
end


mutable struct Segment <: AbstractSegment
    name::String
    id::Int
    residues::Vector{Residue}
    topology::AbstractTopology
    Segment(name::String, id::Int) = begin
        new(name, id, Residue[])
    end
end


# GraphNode{T}(item::T, parent::Opt{GraphNode{T}}) where {T<:AbstractIdentifiable} = begin
#     item.node = GraphNode{T}(item, parent, T[], false, nothing)
# end

const AtomGraphNode = GraphNode{Atom}
const ResidueGraphNode = GraphNode{Residue}


mutable struct Topology <: AbstractTopology
    name::String
    id::Int
    segments::Vector{Segment}
    Topology(name::String, id::Int) = new(name, id, Segment[])
end

# #export IterateByAtom
# struct IterateByAtom
#     target::Topology
#     nseg::Int
#     nres::Int
#     natm::Int
#     IterateByAtom(t::Topology) = begin
#         nseg = length(t.segments)
#         nres = mapreduce(s->length(s.residues), +, t.segments)
#         natm = mapreduce(s->mapreduce(r->length(r.atoms), +, s.residues), +, t.segments)
#         new(t, nseg, nres, natm)
#     end
# end
# Base.iterate(iter::IterateByAtom, (s,r,a)=(1,1,1)) = begin
#     t = iter.target
#     if s > length(t.segments)
#         return nothing
#     elseif r > length(t.segments[s].residues)
#         return iterate(iter, (s+1,1,1))
#     elseif a > length(t.segments[s].residues[r].atoms)
#         return iterate(iter, (s,r+1,1))
#     end
#     (t.segments[s].residues[r].atoms[a], (s,r,a+1))
# end
# Base.size(iter::IterateByAtom) = (iter.nseg, iter.nres, iter.natm)
# Base.length(iter::IterateByAtom) = iter.natm

# Base.length(t::Topology) = begin
#     mapreduce(s->mapreduce(r->length(r.atoms), +, s.residues), +, t.segments)
# end

#---------------------------------------
struct AtomIterator{T <: AbstractIdentifiable}
    target::T
    size::Tuple{Vararg{Int}}
end

AtomIterator(t::Topology) = begin
    nseg = length(t.segments)
    nres = mapreduce(s->length(s.residues), +, t.segments)
    natm = mapreduce(s->mapreduce(r->length(r.atoms), +, s.residues), +, t.segments)
    AtomIterator{Topology}(t, (nseg,nres,natm))
end

AtomIterator(s::Segment) = begin
    nres = length(s.residues)
    natm = mapreduce(r->length(r.atoms), +, s.residues)
    AtomIterator{Segment}(s, (nres,natm))
end

AtomIterator(r::Residue) = begin
    natm = length(r.atoms)
    AtomIterator{Residue}(r, (natm,))
end

Base.size(iter::AtomIterator{T}) where {T<:AbstractIdentifiable} = iter.size
Base.length(iter::AtomIterator{T}) where {T<:AbstractIdentifiable} = iter.size[end]

Base.iterate(iter::AtomIterator{Topology}, (s,r,a)=(1,1,1)) = begin
    t = iter.target
    if s > length(t.segments)
        return nothing
    elseif r > length(t.segments[s].residues)
        return iterate(iter, (s+1,1,1))
    elseif a > length(t.segments[s].residues[r].atoms)
        return iterate(iter, (s,r+1,1))
    end
    (t.segments[s].residues[r].atoms[a], (s,r,a+1))
end

Base.iterate(iter::AtomIterator{Segment}, (r,a)=(1,1)) = begin
    s = iter.target
    if r > length(s.residues)
        return nothing
    elseif a > length(s.residues[r].atoms)
        return iterate(iter, (r+1,1))
    end
    (s.residues[r].atoms[a], (r,a+1))
end

Base.iterate(iter::AtomIterator{Residue}, (a,)=(1,)) = begin
    r = iter.target
    if a > length(r.atoms)
        return nothing
    end
    (r.atoms[a], (a+1,))
end



Base.show(io::IO, node::GraphNode{T}) where {T<:AbstractIdentifiable} = begin
    #println("GraphNode{$(nameof(T))}")
    #println(io, " ├ item:     $(node.item)")
    #println(io, " ├ parent:   $(node.parent===nothing ? "n.d." : node.parent.item)")
    #println(io, " └ children: $(length(node.children))-element array")
    println("GraphNode{$(node.item), $(node.parent===nothing ? "nothing" : node.parent.item)} with $(length(node.children)) children")
end

# Base.reset(r::Residue) = begin
#     foreach(at->at.node.visited=false, r.atoms)
# end

# Base.reset(s::Segment) = begin
#     foreach(r->r.node.visited=false, s.residues)
# end

Base.show(io::IO, at::T) where {T<:AbstractIdentifiable} = begin
    print(io, "$(nameof(T)){$(at.name):$(at.id)}")
end

Base.push!(top::Topology, seg::Segment) = begin
    push!(top.segments, seg)
    seg.topology = top
    top
end
Base.push!(seg::Segment, res::Residue) = begin
    push!(seg.residues, res)
    res.segment = seg
    seg
end
Base.push!(res::Residue, atom::Atom) = begin
    res.atomsbyname[atom.name] = atom
    push!(res.atoms, atom)
    atom.residue = res
    res
end

Base.get(res::Residue, name::String, default) = get(res.atomsbyname, name, default)
Base.get(res::Residue, name::String) = get(res, name, nothing)

Base.push!(parent::GraphNode{T}, node::GraphNode{T}) where {T<:AbstractIdentifiable} = begin
    node.parent = parent
    push!(parent.children, node)
end

export  isorphan
isorphan(node::GraphNode{T}) where {T<:AbstractIdentifiable} = begin
    node.parent===nothing && isempty(node.children)
end

# Base.split(r1::GraphNode{T}, r2::GraphNode{T}) where {T<:AbstractIdentifiable} = begin
#     if r2.parent === r1
#         i = findfirst(r->r===r2, r1.children)
#         i !== nothing && deleteat!(r1.children, i)
#         r2.parent = nothing
#     elseif r1.parent === r2
#         split(r2, r1)
#     end
# end

export hasparent
# hasparent(n::GraphNode{T}) where {T<:AbstractIdentifiable} = n.parent !== nothing
hasparent(n::GraphNode) = n.parent !== nothing
hasparent(::Nothing) = false


Base.show(io::IO, at::Atom) = begin
    if isdefined(at, :residue) && 
        isdefined(at.residue, :segment) &&
        isdefined(at.residue.segment, :topology)
        r = at.residue
        s = r.segment
        t = s.topology
        print(io, "Atom/$(t.name):$(t.id)/$(s.name):$(s.id)/$(r.name):$(r.id)/$(at.name):$(at.id)")
    elseif isdefined(at, :residue) && 
        isdefined(at.residue, :segment)
        r = at.residue
        s = r.segment
        print(io, "Atom/undef/$(s.name):$(s.id)/$(r.name):$(r.id)/$(at.name):$(at.id)")
    elseif isdefined(at, :residue)
        r = at.residue
        print(io, "Atom/undef/undef/$(r.name):$(r.id)/$(at.name):$(at.id)")
    else
        print(io, "Atom/undef/undef/undef/$(at.name):$(at.id)")
    end
end

# export iterbysegment, iterbyresidue, iterbyatom
# # iterbyatom(t::Topology) = (a for a in t.atoms)
# iterbyatom(t::Topology) = IterateByAtom(t)
# iterbyatom(r::Residue) = (a for a in r.atoms)
# iterbyresidue(t::Topology) = (r for s in t.segments for r in s.residues)
# iterbysegment(t::Topology) = (s for s in t.segments)

export eachatom
eachatom(t::T) where {T <: AbstractIdentifiable} = AtomIterator(t)

export traverse
traverse(f::F, node::GraphNode{T}) where {T<:AbstractIdentifiable, F<:Function} = begin
    queue = Vector{GraphNode{T}}()
    push!(queue, node)
    while !isempty(queue)
        pivot = popfirst!(queue)
        for node in pivot.children
            push!(queue, node)
        end
        f(pivot)
    end
end

export ref
const ref = Residue("ref", 0)
let
    at1 = Atom("Y", -1, -2, "?")
    at2 = Atom("X", -2, -1, "?")
    at3 = Atom("O", -3,  0, "?")
    push!(at1.node, at2.node)
    push!(at2.node, at3.node)
    push!(ref, at1)
    push!(ref, at2)
    push!(ref, at3)
end

using LinearAlgebra: I

export NodeState
mutable struct NodeState{T<:AbstractFloat}
    # t::Vector{Float64}  # translation vector
    # r::Matrix{Float64}  # local to global rotation matrix
    
    t::MVector{3,T}  # translation vector
    r::MMatrix{3,3,T,9}  # local to global rotation matrix
    
    # internal coordinates
    b::T  # self<->parent bond length
    θ::T  # self<->parent<->grandparent angle
    ϕ::T  # self<->parent<->grandparent<->grand-grandparent dihedral
    
    Δϕ::T
    changed::Bool
    #NodeState{T}() where {T<:AbstractFloat} = begin
    #    new(zeros(3), Matrix(1.0I,3,3), 0.0, 0.0, 0.0)
    #end
end
NodeState{T}() where {T} = begin
    NodeState{T}(
        MVector{3,T}(zeros(T, 3)),
        MMatrix{3,3,T,9}(Matrix{T}(I,3,3)),
        zero(T),
        zero(T),
        zero(T),
        zero(T),
        false
    )
end

Base.setproperty!(ns::NodeState{T}, key::Symbol, val) where T = begin
    if key == :Δϕ
        setfield!(ns, :changed, true)
    end
    setfield!(ns, key, val)
end

export State
mutable struct State{T<:AbstractFloat}
    items::Vector{NodeState{T}}
    size::Int
    id::Int
    i2c::Bool
    c2i::Bool
    x::Array{T,2}
    # f::Array{T,2}
    # v::Array{T,2}
end
State{T}(n::Int) where T = begin
    items = Vector{NodeState{T}}(undef, n+3)
    for i=1:n+3
        items[i] = NodeState{T}()
    end
    items[1].t[1] = -1.0
    items[1].t[2] =  1.0
    items[2].t[1] = -1.0
    State{T}(items, n, -1, false, false, zeros(T,3,n))
end
State(::Type{T}, n::Int) where T = State{T}(n)
State(n::Int) = State{Float64}(n)

Base.getindex(s::State, i::Int) = begin
    s.items[i+3]
end

Base.firstindex(s::State) = -2
Base.lastindex(s::State) = s.size

Base.eltype(::Type{State{T}}) where T = T







const PDB = Val{1}


Base.read(filename::AbstractString, ::Type{PDB}) = begin
    top,state = open(filename) do fin
        read(fin, PDB)
    end
    top.name = basename(filename)
    top,state
end

Base.read(io::IO, ::Type{PDB}) = begin
    
    top = Topology("?", -1)
    seg =  Segment("?", -1)     # orphan segment
    res =  Residue("?", -1)     # orphan residue
    
    seekstart(io)
    natoms = mapreduce(l->startswith(l, "ATOM")||startswith(l, "HETATM"), +, eachline(io); init=0)
    
    state = State(natoms)
    # request conversion from cartesian to internal
    state.c2i = true

    id2atom = Dict{Int,Atom}()
    segid = atmindex = 1
    
    seekstart(io)
    for line in eachline(io)
        
        if startswith(line, "ATOM") || startswith(line, "HETATM")
            
            resname = string(strip(line[18:20]))
            segname = string(strip(string(line[22])))
            resid = parse(Int, line[23:26])

            if seg.name != segname
                seg = Segment(segname, segid)
                push!(top, seg)
                segid += 1
            end

            if res.id != resid || res.name != resname
                res = Residue(resname, resid)
                push!(seg, res)
            end

            atsymbol = length(line)>77 ? string(strip(line[77:78])) : "?"
            atname = string(strip(line[13:16]))
            atid = parse(Int, line[7:11])

            atom = Atom(atname, atid, atmindex, atsymbol)
            id2atom[atid] = atom
            push!(res, atom)
            
            s = state[atmindex]
            s.t[1] = parse(Float64, line[31:38])
            s.t[2] = parse(Float64, line[39:46])
            s.t[3] = parse(Float64, line[47:54])
            atmindex += 1

        elseif startswith(line, "CONECT")
            idxs = map(s->parse(Int, s), [line[n:n+4] for n=7:5:length(line)])
            pivot = id2atom[idxs[1]]
            for i in idxs[2:end]
                other = id2atom[i]
                other ∉ pivot.bonds && push!(pivot.bonds, other)
                pivot ∉ other.bonds && push!(other.bonds, pivot)
            end
        end
    end
    top, state
end

using Printf

# Base.write(io::IO, top::Topology, state::State, ::Type{Zmat.PDB}) = begin
Base.write(io::IO, top::Topology, state::State) = begin
    
    println(io, "MODEL")
    for atom in eachatom(top)
        istate = state[atom.index]
        s = @sprintf("ATOM %6d %4s %-4sA %3d    %8.3f%8.3f%8.3f%24s",
            atom.index, atom.name,
            atom.residue.name, atom.residue.id,
            istate.t[1],
            istate.t[2],
            istate.t[3],
            atom.symbol)
        println(io, s)
    end

    for atom in eachatom(top)
        print(io, @sprintf("CONECT%5d", atom.index))
        foreach(n->print(io, @sprintf("%5d",n.item.index)), atom.node.children)
        println(io,"")
    end
    println(io, "ENDMDL")
end

export build_tree!
function build_tree!(top::Topology)
    queue = Atom[]
    byatom = eachatom(top)
    
    # reset graph flags
    foreach(at->at.node.visited=false, byatom)
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

    root = get(ref, "O").node
    for atom in byatom
        node = atom.node

        # if this atom is orphan, then make it a child
        # of the root (origin)
        if !hasparent(node)
            push!(root, node)
           continue
        end

        # build residue graph
        child = node.item.residue
        parent = node.parent.item.residue
        if child===parent || (child.node.visited && parent.node.visited)
            continue
        end
        push!(parent.node, child.node)
        parent.node.visited = true
        child.node.visited = true
    end

    for atom in byatom
        node = atom.node
        node.ascendents = (
            node.item.index,
            node.parent.item.index,
            node.parent.parent.item.index,
            node.parent.parent.parent.item.index
        )
    end

end


export sync!
sync!(state::State, top::Topology) = begin
    if state.c2i && state.i2c
        throw(Exception)
    elseif state.c2i
        c2i!(state, top)
    elseif state.i2c
        i2c!(state, top)
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
        # i = atom.index
        # j = atom.node.parent.item.index
        # k = atom.node.parent.parent.item.index
        # l = atom.node.parent.parent.parent.item.index
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
    
    queue = AtomGraphNode[]

    xyz = zeros(3, state.size)
    root = get(ref, "O").node
    append!(queue, root.children)
    
    while !isempty(queue)
        node = popfirst!(queue)
        (i,j,k) = node.ascendents
        # i = node.item.index
        istate = state[i]
        
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
        Zmat.@cross u n[u] vji[u] vjk[u]
        dn = sqrt(dot(n,n))
        @nexprs 3 u -> Ri[u,3] = n[u]/dn
    
        # column 2 (y)
        Zmat.@cross u Ri[u,2] Ri[u,3] Ri[u,1]
        
        # move to new position
        # @nexprs 3 u -> istate.t[u] = vji_u + jstate.t[u]
        # xyz[:, i] .= istate.t
        @. istate.t = vji + jstate.t
        @. xyz[:,i] = istate.t
    end
    xyz
end



end