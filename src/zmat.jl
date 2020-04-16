module Treta

include("ProtoSyn.jl")

ala = Dict(
    "N"   => ("-C", "-CA",  "-N"),
    "H"   => ( "N",  "-C", "-CA"),
    "CA"  => ( "N",  "-C", "-CA"),
    "C"   => ("CA",   "N",  "-C"),
    "O"   => ( "C",  "CA",   "N"),
    "HA"  => ("CA",   "N",  "-C"),
    "CB"  => ("CA",   "N",  "-C"),
    "HB1" => ("CB",  "CA",   "N"),
    "HB2" => ("CB",  "CA",   "N"),
    "HB3" => ("CB",  "CA",   "N")
)

lib = ProtoSyn.Peptides.load(ProtoSyn.Peptides.resource_dir)

mol, state = ProtoSyn.Peptides.build("FF", lib)

pymol = ProtoSyn.XMLRPC.ServerProxy("http://localhost", 9123)
pymol.delete("all")

function render(mol::ProtoSyn.Molecule, state::ProtoSyn.State, name::String)
    io = IOBuffer()
    write(io, mol, state, ProtoSyn.PDB)
    pymol.read_pdbstr(String(take!(io)), name)
end

render(mol, state, mol.name)

function traverse2(pivot::Int, path::Tuple, depth::Int, visited::BitArray, bonds::ProtoSyn.BondGraph)
    #println("x0:", depth, path, pivot)
    path = tuple(path..., pivot)
    #println("x1:", depth, path)
    if depth == 3
        println("=> ", path)
    else
        visited[pivot] = true
        for i in bonds[pivot]
            if !visited[i]
                traverse2(i, path, depth+1, visited, bonds)
            end
        end
        visited[pivot] = false
    end
end

vis = falses(size(mol))
# p = get(mol.residues[1], "N")
# println(p)
# traverse2(p.id, (), 0, vis, mol.bonds)

for atom in ProtoSyn.iterbyatom(mol)
    id = atom.id + atom.residue.offset
    if length(mol.bonds[id]) > 1
        traverse2(id, (), 0, vis, mol.bonds)
    end
    vis[id] = true
end


struct TreeNode
    atom::ProtoSyn.Atom
    parent::ProtoSyn.Opt{TreeNode}
    children::Vector{TreeNode}
end

TreeNode(at::ProtoSyn.Atom, parent::ProtoSyn.Opt{TreeNode}) = TreeNode(at, parent, [])
TreeNode(at::ProtoSyn.Atom) = TreeNode(at, nothing, [])

Base.show(io::IO, node::TreeNode) = begin
    at = node.atom
    children = map(n->n.atom.id+n.atom.residue.offset, node.children)
    println(io, "TreeNode")
    println(io, " ├ atom: /$(at.residue.molecule.name)/$(at.residue.name)/$(at.name)")
    println(io, " ├ parent: $(node.parent.atom.id+node.parent.atom.residue.offset)")
    println(io, " └ children: $children")    
end


struct Topology
    segments::Vector{TreeNode}
    atoms::Vector{TreeNode}
end

struct Residue
    name::String
    id::Int
end



struct _Object{T}
    name::String
    id::Int
    chains::Dict{String, T}
end

struct _Segment{T}
    name::String
    id::Int
    object::_Object
    residues::Vector{T}
end

struct _Residue{T}
    name::String
    id::Int
    segment::_Segment
    atoms::Vector{T}
end

struct Atom
    name::String
    id::Int
    index::Int
    # resname::String
    # resid::Int
    # segname::String
    # segid::Int
    residue::_Residue
end

const Residue = _Residue{Atom}
const Segment = _Segment{Residue}
const Object = _Object{Segment}

Base.show(io::IO, at::Atom) = begin
    res = at.residue
    seg = res.segment
    obj = seg.object
    # print(io, "Atom{/$(at.segname):$(at.segid)/$(at.resname):$(at.resid)/$(at.name):$(at.id)}")    
    print(io, "Atom{/$(obj.name):$(obj.id)/$(seg.name):$(seg.id)/$(res.name):$(res.id)/$(at.name):$(at.id)}")
end

obj = Object("obj", 1, Dict())
seg = Segment("A", 1, obj, [])
obj.chains[seg.name] = seg

let res = nothing
for (index,atom) in enumerate(ProtoSyn.iterbyatom(mol))
    if isa(res, Nothing) || (atom.residue.name, atom.residue.id) != (res.name, res.id)
        res = Residue(atom.residue.name, atom.residue.id, seg, [])
        push!(seg.residues, res)
    end
    at = Atom(atom.name, atom.id, index, res)
    push!(res.atoms, at)
end
end

println(obj)

# atoms = map(
#     at->Atom(at.name,at.id,at.id+at.residue.offset,at.residue.name,at.residue.id,"",1),
#     ProtoSyn.iterbyatom(mol)
# )

# println(atoms)

#mutable 


function build_tree(molecule::ProtoSyn.Molecule)
    ts = ProtoSyn.TraverseState(size(molecule))
    head = tail = 0

    root = TreeNode(get(molecule.residues[1], "N"))
    tree = Vector{TreeNode}(undef, size(molecule))
    
    ts.visited[root.atom.id+root.atom.residue.offset] = true
    tree[tail+=1] = root

    #root.id+root.residue.offset
    
    atoms = collect(ProtoSyn.iterbyatom(molecule))

    while head < tail
        # id1 = ts.indices[head+=1]
        parent = tree[head+=1]
        id1 = parent.atom.id + parent.atom.residue.offset
        for id2 in molecule.bonds[id1]
            if !ts.visited[id2]
                node = TreeNode(atoms[id2], root)
                push!(parent.children, node)
                ts.visited[id2] = true
                # ts.indices[tail+=1] = id2
                tree[tail+=1] = node
            end
        end
    end
    return tree
end

tree = build_tree(mol)

end