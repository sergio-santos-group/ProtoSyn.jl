module Treta

include("src/ProtoSyn.jl")


lib = ProtoSyn.Peptides.load(ProtoSyn.Peptides.resource_dir)

mol, state = ProtoSyn.Peptides.build("F", lib)
mol2, state2 = ProtoSyn.Peptides.build("AGP", lib)

pymol = ProtoSyn.XMLRPC.ServerProxy("http://localhost", 9123)
# pymol.delete("all")

function render(mol::ProtoSyn.Molecule, state::ProtoSyn.State, name::String)
    io = IOBuffer()
    write(io, mol, state, ProtoSyn.PDB)
    pymol.read_pdbstr(String(take!(io)), name)
end

# render(mol, state, mol.name)




abstract type AbstractNode end
abstract type AbstractTopology end

mutable struct _Atom{T<:AbstractNode}
    name::String
    id::Int
    index::Int
    symbol::String
    resname::String
    resid::Int
    segname::String
    segid::Int
    bonds::Vector{_Atom}
    node::T
    function _Atom{T}(name::String,
        id::Int,
        index::Int,
        symbol::String,
        resname::String,
        resid::Int,
        segname::String,
        segid::Int) where {T<:AbstractNode}
        return new{T}(name,id,index,symbol,resname,resid,segname,segid,[])
    end
end


mutable struct TreeNode <: AbstractNode
    atom::_Atom
    parent::Union{Nothing,TreeNode}
    children::Vector{TreeNode}
    # rmat::Matrix{Float64}
    visited::Bool
    topology::AbstractTopology
    TreeNode(at::_Atom, p::Union{Nothing,TreeNode}) = begin
        node = new(at,p,[],false)
        at.node = node
        return node
    end
end
const Atom = _Atom{TreeNode}

Base.show(io::IO, at::Atom) = begin
    print(io, "Atom/$(at.segname)-$(at.segid)/$(at.resname)-$(at.resid)/$(at.name)-$(at.id)#$(at.index)")    
end

#TreeNode(at::Atom, parent::TreeNode) = TreeNode(at, parent)
TreeNode(at::Atom) = TreeNode(at, nothing)

Base.show(io::IO, node::TreeNode) = begin
    println(io, "TreeNode")
    println(io, " ├ atom: $(node.atom)")
    println(io, " ├ parent: $(node.parent===nothing ? "nothing" : node.parent.atom)")
    println(io, " └ children: $(length(node.children))-element array")
end

mutable struct Residue
    name::String
    id::Int
    atoms::Atom
    parent::
end

mutable struct Segment
    name::String
    id::Int
    residues::Residue
end

# mutable struct Topology{T,K}
mutable struct Topology <: AbstractTopology
    nodes::Vector{TreeNode}
    atoms::Vector{Atom}
    
    _indexbyresid::Dict{Int,Residue}
    _indexbysegid::Dict{Int,Vector{Atom}}
    
    index_offset::Int
    resid_offset::Int
    size::Int
end



# Topology(size::Int) = Topology(Vector{TreeNode}(undef, size), Vector{Atom}(undef, size), size, Dict(), Dict())
Topology() = Topology([], [], Dict(), Dict(), 0, 0, 0)
# Topology(size::Int) = Topology([], Vector{TreeNode}(undef, size), Vector{Atom}(undef, size), size, Dict())

Base.iterate(t::Topology, idx=1) = idx > t.size ? nothing : (t.atoms[idx], idx+1)
Base.length(t::Topology) = t.size
# Base.size(t::Topology) = t.size

Base.push!(t::Topology, n::TreeNode) = begin
    push!(t.nodes, n)
    n.topology = t
    return t
end

Base.push!(t::Topology, at::Atom) = begin
    at.index = (t.index_offset += 1)
    push!(t.atoms, at)
    t.size +=1 ;
    ar = get!(t._indexbyresid, at.resid) do; []; end
    push!(ar, at)
    ar = get!(t._indexbysegid, at.segid) do; []; end
    push!(ar, at)
    return t
end

# top = Topology([],[],[])
# top = Topology()
# top2 = Topology()

# top = Topology(size(mol))

function mol2top(m::ProtoSyn.Molecule)
    top = Topology()
    for at in ProtoSyn.iterbyatom(m)
        atom = Atom(at.name, at.id, 0, at.symbol,
                at.residue.name, at.residue.id,
                "", 1)
        push!(top, atom)
    end
    for (pivot, others) in m.bonds
        top.atoms[pivot].bonds = top.atoms[others]
    end
    top
end

# # ADD ATOMS
# for (index,at) in enumerate(ProtoSyn.iterbyatom(mol))
#     atom = Atom(at.name, at.id, 0, at.symbol,
#                 at.residue.name, at.residue.id,
#                 "", 1)
#     push!(top, atom)
#     atom = Atom(at.name, at.id, 0, at.symbol,
#                 at.residue.name, at.residue.id,
#                 "", 1)
#     push!(top2, atom)
# end
# # ADD BONDS
# for (pivot, others) in mol.bonds
#     top.atoms[pivot].bonds = top.atoms[others]
#     top2.atoms[pivot].bonds = top2.atoms[others]
# end

struct Selection
    segid::Int
    resid::Int
    atom::AbstractString
end
Selection(str::String) = begin
    fields = split(str, "/")
    Selection(parse(Int,fields[2]), parse(Int,fields[3]), fields[4])
end

function (s::Selection)(t::Topology)
    get(t,s,nothing)
end

macro s_str(str)
    Selection(str)
end

Base.get(t::Topology, s::Selection, default) = begin
    container = get(t._indexbyresid, s.resid, nothing)
    if container===nothing
        return default
    end
    i = findfirst(at->at.name==s.atom, container)
    if i===nothing
        return default
    end
    i===nothing ? default : container[i]
end


Base.merge!(combine::Function, res1::Residue, res2::Residue) = begin
    
    top1 = res1[1].node.topology
    top2 = res2[1].node.topology
    
    # is this an intra-topology merge?
    if top1===top2
        combine(res1, res2)
    else
        top2 = deepcopy(top2)
        com
    end

    other = top1!==top2 ? deepcopy(top2) : top2
    combine(top1, other)


    for atom in other.atoms
        push!(dst, atom)
        push!(dst, atom.node)
    end
    dst
end
# Base.merge!(combine::Function, dst::Topology, other::Topology) = begin
#     other = deepcopy(other)
#     combine(dst, other)
#     for atom in other.atoms
#         push!(dst, atom)
#         push!(dst, atom.node)
#     end
#     dst
# end

function build_tree!(topology::Topology)
    
    visited = falses(topology.size)
    head = tail = 0

    for i=1:topology.size

        visited[i] && continue

        root = TreeNode(topology.atoms[i])
        
        visited[root.atom.id] = true
        push!(topology, root)
        tail += 1

        while head < tail
            parent = topology.nodes[head+=1]
            for atom in parent.atom.bonds
                if !visited[atom.index]
                    node = TreeNode(atom, parent)
                    push!(parent.children, node)
                    visited[atom.index] = true
                    push!(topology, node)
                    tail += 1
                end
            end
        end
    end # end for

    return topology
end

function connect(at1::Atom, at2::Atom)
    at1 ∉ at2.bonds && push!(at2.bonds, at1)
    at2 ∉ at1.bonds && push!(at1.bonds, at2)
    parent = at1.node
    child = at2.node
    child ∉ parent.children && push!(parent.children, child)
    child.parent = parent
end

# PEPTIDE SPECIFIC MERGE
function peptide_merge(t1::Topology, t2::Topology)

    C = s"/1/1/C"(t1)
    N = s"/1/1/N"(t2)
    
    for node in t2.nodes
        node.visited = false
    end
    
    # update all residue names and ids
    # in this molecule (all atoms reachable
    # from atom N)
    queue = TreeNode[N.node]
    N.node.visited = true
    d = Dict{Int,Int}()
    resid = max(
        maximum(keys(t1._indexbyresid)),
        maximum(keys(t2._indexbyresid))
    )

    while !isempty(queue)
        node = popfirst!(queue)

        node.atom.segid = C.segid
        node.atom.segname = C.segname
        node.atom.resid = get!(d,node.atom.resid) do; resid+=1; end        

        for other in node.atom.bonds
            if !other.node.visited
                other.node.visited = true
                push!(queue, other.node)
            end
        end
    end
    
    connect(C, N)
    
end

top = mol2top(mol)
top2 = mol2top(mol2)
build_tree!(top)
build_tree!(top2)

merge!(peptide_merge, top, top2)


end