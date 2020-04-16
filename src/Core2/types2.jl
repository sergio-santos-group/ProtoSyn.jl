#module TT

# optional type
export Opt
const Opt = Union{Nothing, T} where T


export Atom, Residue, Molecule, Link
#export Link, LinkedResidue
#export State, Dihedral, AxisRotatableBlock


export ResidueLib

# connectivity graph
export BondGraph
const BondGraph = Dict{Int, Vector{Int}}
#const BondGraphByName = Dict{String, Vector{String}}




export TraverseState
mutable struct TraverseState
    visited::BitVector
    indices::Vector{Int}
    size::Int
    capacity::Int
    TraverseState(n::Int) = begin
        new(falses(n), zeros(Int, n), 0, n)
    end
end


# ===================================================================
mutable struct _Link{T}
    r1::T
    r2::T
    bonds::BondGraph
    flagged::Bool
end

# ===================================================================
mutable struct _Molecule{T}
    id::Int
    name::String
    residues::Vector{T}
    root::Opt{T}
    bonds::BondGraph
    links::Vector{_Link{T}}
    valid::Bool
    traverse_state::Opt{TraverseState}
end



# ===================================================================
mutable struct _Residue{T}
    id::Int
    name::String
    molecule::Opt{_Molecule}
    # chain::Opt{_Chain}
    atoms::Vector{T}
    atomsbyname::Dict{String,T}
    bonds::BondGraph
    valid::Bool
    links::Vector{_Link{_Residue{T}}}
    offset::Int
    flagged::Bool
    #metadata::Opt{Dict}

    coords::Opt{AbstractArray}
    # coords::SubArray
end

# ===================================================================
mutable struct Atom
    id::Int
    #lid::Int
    name::String
    symbol::String
    residue::Opt{_Residue}
    metadata::Opt{Dict}
end



const Residue = _Residue{Atom}
const Molecule = _Molecule{Residue}
const Link = _Link{Residue}
const ResidueLib = Dict{String, Residue}

Link(r1::Residue, r2::Residue) = Link(r1,r2,BondGraph(),false)
Link(r1::Residue, r2::Residue, pairs::Pair{Int,Int}...) = begin
    l = Link(r1, r2, BondGraph(), false)
    foreach(p->push!(l,p), pairs)
    l
end

Atom() = Atom(0, "UNKATM", "?", nothing, nothing)
Residue() = Residue(0, "UNKRES", nothing, Atom[], Dict{String,Atom}(), BondGraph(), true, Link[], 0, false, nothing)
Molecule() = Molecule(0, "UNKMOL", Residue[], nothing, BondGraph(), Link[], true, nothing)
# Residue(m::Molecule) = begin
#     r = Residue()
#     push!(m, r)
#     r
# end

# Atom(r::Residue) = begin
#     a = Atom()
#     push!(r, a)
#     a
# end
# ===================================================================
export AbstractAxisRotatableBlock, Dihedral, AxisRotatableBlock

abstract type AbstractAxisRotatableBlock end

struct Dihedral <: AbstractAxisRotatableBlock
    a0::Int
    a1::Int
    a2::Int
    a3::Int
    root::Int
    type::Opt{Enum{<:Integer}}
end
Dihedral(a0::Int, a1::Int, a2::Int, a3::Int) = Dihedral(a0, a1, a2, a3, a2, nothing)
Dihedral(a0::Int, a1::Int, a2::Int, a3::Int, type::Opt{Enum{<:Integer}}) = 
    Dihedral(a0, a1, a2, a3, a2, type)

Dihedral(a0::Atom, a1::Atom, a2::Atom, a3::Atom, type::Opt{Enum{<:Integer}}=nothing) = 
    Dihedral(
        a0.id+a0.residue.offset,
        a1.id+a1.residue.offset,
        a2.id+a2.residue.offset,
        a3.id+a3.residue.offset,
        a2.id+a2.residue.offset,
        type)

#@kwdef
struct AxisRotatableBlock <: AbstractAxisRotatableBlock
    a1::Int
    a2::Int
    root::Int
    #enpoints::Vector{Int} = []
    # movable::Opt{Vector{Int}} = nothing
    # terminals::Opt{Vector{Int}} = nothing
end



# ===================================================================


Base.isvalid(mol::Molecule) = mol.valid

invalidate(res::Residue) = begin
    if res.valid
        res.valid = false
        res.molecule !== nothing && invalidate(res.molecule)
    end
end


invalidate(mol::Molecule) = begin
    mol.valid = false
end


function update!(mol::Molecule)

    if isempty(mol.residues)
        mol.valid = true
        return mol
    end

    if mol.root === nothing
        error("molecule update requires specification of a root residue")
    end

    # # sort residues by breadth-first traversal 
    # foreach(r -> r.flagged=false, mol.residues)
    # head = 0
    # tail = 1
    # mol.root.flagged = true
    # mol.residues[tail] = mol.root
    # while head < tail
    #     res = mol.residues[head+=1]
    #     for l in res.links
    #         other = l.r1===res ? l.r2 : l.r1
    #         if !other.flagged
    #             mol.residues[tail+=1] = other
    #             other.flagged = true
    #         end
    #     end
    # end

    # if any(r->!r.flagged, mol.residues)
    #     error("free floating residues found within molecule")
    # end

    # # renumber residues and atoms
    # atmid = resid = 0
    # for res in mol.residues
    #     #lid = 0
    #     res.id = (resid += 1)
    #     res.offset = atmid
    #     for atom in res.atoms
    #         atom.id = (atmid += 1)
    #         atom.lid = (lid += 1)
    #     end
    # end

    # regenerate bond graph (at molecular level)
    #  1st pass deals only with intra-residue bonds
    offset = resid = 0
    for res in mol.residues
        #res.id = (resid += 1)
        res.offset = offset
        for (pivot,others) in res.bonds
            mol.bonds[pivot+offset] = offset .+ others
        end
        offset += length(res)
    end

    # second pass now deals with inter-residue links
    for link in mol.links
        offset1 = link.r1.offset
        offset2 = link.r2.offset
        for atm1 in keys(link.bonds)
            for atm2 in link.bonds[atm1]
                push!(get!(mol.bonds, atm1+offset1, []), atm2+offset2)
                push!(get!(mol.bonds, atm2+offset2, []), atm1+offset1)
            end
        end

    end

    natoms = size(mol)
    if (mol.traverse_state===nothing) || (natoms > mol.traverse_state.capacity)
        mol.traverse_state = TraverseState(natoms)
    end

    mol.valid = true
    return mol

end



function update!(ts::TraverseState, rotblock::AbstractAxisRotatableBlock, bonds::BondGraph)
    head = tail = 0
    
    fill!(ts.visited, false)
    #ts.indices[tail+=1] = rotblock.a1
    ts.visited[rotblock.a1] = true
    ts.visited[rotblock.a2] = true

    ts.indices[tail+=1] = rotblock.root
    while head < tail
        id1 = ts.indices[head+=1]
        for id2 in bonds[id1]
            if !ts.visited[id2]
                ts.visited[id2] = true
                ts.indices[tail+=1] = id2
            end
        end
    end
    ts.size = tail
    return ts
end

# @inline Base.reset(ts::TraverseState) = begin
#     fill!(ts.visited, false)
#     ts
# end

Base.push!(residue::Residue, atom::Atom) = begin

    if atom.residue!==nothing
        error("only orphan atoms can be added to a residue.")
    end
    if haskey(residue.atomsbyname, atom.name)
        error("Atom names within residues must be unique.")
    end

    atom.residue = residue
    push!(residue.atoms, atom)
    residue.atomsbyname[atom.name] = atom

    invalidate(residue)
    residue
end


Base.push!(molecule::Molecule, residue::Residue) = begin

    if residue.molecule!==nothing
        error("only orphan residues can be added to a molecule.")
    end
    if !in(residue, molecule.residues)
        invalidate(molecule)
        residue.molecule = molecule
        push!(molecule.residues, residue)
    end
    molecule
end

# Base.delete!(residue::Residue, link::Link) = begin
#     idx = findfirst(l->l===link, residue.links)
#     if idx !== nothing
#         deleteat!(residue.links, idx)
#     end
#     residue
# end

Base.delete!(container::Union{Residue,Molecule}, link::Link) = begin
    idx = findfirst(l->l===link, container.links)
    if idx !== nothing
        deleteat!(container.links, idx)
    end
    container
end

Base.delete!(molecule::Molecule, residue::Residue) = begin
    if residue.molecule !== molecule
        error("unable to delete residues belonging to other molecules")
    end
    idx = findfirst(r->r===residue, molecule.residues)
    if idx === nothing
        return molecule
    end

    deleteat!(molecule.residues, idx)
    residue.molecule = nothing
    invalidate(molecule)
    
    # remove links
    for link in residue.links
        # remove from molecule
        delete!(molecule, link)
        
        # remove link from adjacent residues
        other = link.r1!==residue ? link.r1 : link.r2
        delete!(other, link)
    end
    molecule
end
Base.pop!(molecule::Molecule, residue::Residue) = begin
    delete!(molecule, residue)
    residue
end

Base.push!(l::Link, pair::Pair{Int,Int}) = begin
    if !haskey(l.bonds, pair.first)
        l.bonds[pair.first] = [pair.second]
    else
        push!(l.bonds[pair.first], pair.second)
    end
    l
end

Base.in(r::Residue, m::Molecule) = (r.molecule === m) || Base.in(r, m.residues)
Base.in(a::Atom, r::Residue) = (a.residue === r) || Base.in(a, r.atoms)
Base.in(s::AbstractString, r::Residue) = haskey(r.atomsbyname, s)


Base.get(r::Residue, k::AbstractString, default) = get(r.atomsbyname, k, default)
Base.get(r::Residue, k::AbstractString) = get(r.atomsbyname, k, nothing)

Base.length(r::Residue) = length(r.atoms)
Base.length(m::Molecule) = length(m.residues)
Base.size(m::Molecule) = mapreduce(r->length(r), +, m.residues; init=0)
Base.size(r::Residue) = length(r)

# Base.insert!(molecule::Molecule, index::Integer, residue::Residue) = begin
#     if residue.molecule!==nothing
#         error("only orphan residues can be inserted into a molecule.")
#     end
#     insert!(molecule.residues, index, residue)
#     invalidate(molecule)
# end

function link(linkprovider::F, r1::Residue, r2::Residue) where {F<:Function}
    if r1.molecule !== r2.molecule
        error("unable to link residues belonging to different molecules")
    end
    l::Link = linkprovider(r1,r2)
    if l!==nothing
        push!(r1.links, l)
        push!(r2.links, l)
        push!(r1.molecule.links, l)
    end
    l
end

Base.delete!(container::Union{Residue,Molecule}, l::Link) = begin
    if  (idx = findfirst(item->item===l, container.links)) !== nothing
        deleteat!(container.links, idx)
    end
    container
end


function repr_vec(item::AbstractArray)
    "$(size(item,1))-element array"
end

function repr_opt(item)
    "$(item.name).$(item.id)"
end

function repr_opt(item::Nothing)
    "n.d."
end

Base.show(io::IO, item::Atom) = begin
    println(io, "Atom")
    println(io, " ├ id: $(item.id)")
    println(io, " ├ id: $(item.id+item.residue.offset)")
    println(io, " ├ name: $(item.name)")
    println(io, " ├ symbol: $(item.symbol)")
    println(io, " └ residue: $(repr_opt(item.residue))")
end

Base.show(io::IO, item::Residue) = begin
    println(io, "Residue")
    println(io, " ├ id: $(item.id)")
    println(io, " ├ name: $(item.name)")
    println(io, " ├ offset: $(item.offset)")
    println(io, " ├ molecule: $(repr_opt(item.molecule))")
    println(io, " ├ links: $(repr_vec(item.links))")
    println(io, " └ atoms: $(repr_vec(item.atoms))")
end

Base.show(io::IO, item::Molecule) = begin
    println(io, "Molecule")
    println(io, " ├ id: $(item.id)")
    println(io, " ├ name: $(item.name)")
    println(io, " ├ valid: $(item.valid)")
    println(io, " ├ root: $(repr_opt(item.root))")
    println(io, " └ residues: $(repr_vec(item.residues))")
end

Base.show(io::IO, item::Link) = begin
    println(io, "Link")
    println(io, " ├ r1: $(item.r1.name).$(item.r1.id)")
    println(io, " ├ r2: $(item.r2.name).$(item.r2.id)")
    println(io, " └ bonds:")
    for pair in pairs(item.bonds)
        println(io, "   ├ $(pair.first) => $(join(map(string, pair.second), ", "))")
    end
end

export iterbyatom
iterbyatom(r::Residue)  = (a for a in r.atoms)
iterbyatom(m::Molecule) = (a for r in m.residues for a in r.atoms)

export eachatom
eachatom(r::Residue)  = (a for a in r.atoms)
eachatom(m::Molecule) = (a for r in m.residues for a in r.atoms)



istemplate(r::Residue) = begin
    # (r.molecule === nothing && isa(r.metadata, Dict) && haskey(r.metadata,:coords))
    r.molecule === nothing
end


# res4(n) = begin
#     r = Residue()
#     r.name = "square$(n)"
#     r.bonds[1] = [2,3]
#     r.bonds[2] = [1,4]
#     r.bonds[3] = [1,4]
#     r.bonds[4] = [2,3]
#     r.atoms = [Atom() for i in 1:4]
#     return r
# end

# res3(n) = begin
#     r = Residue()
#     r.name = "triangle$(n)"
#     r.bonds[1] = [2,3]
#     r.bonds[2] = [1,3]
#     r.bonds[3] = [1,2]
#     r.atoms = [Atom() for i in 1:3]
#     return r
# end

# r1 = res4(1)
# r2 = res3(2)
# r3 = res3(3)
# r4 = res4(4)
# r5 = res4(5)
# r6 = res4(6)
# mol = Molecule()
# push!(mol, r1)
# push!(mol, r2)
# push!(mol, r3)
# push!(mol, r4)
# push!(mol, r5)
# push!(mol, r6)
# mol.root = r1

# link((a,b)->Link(a,b,BondGraph(3=>[1]), false), r1, r2)
# link((a,b)->Link(a,b,BondGraph(2=>[1]), false), r2, r3)
# link((a,b)->Link(a,b,BondGraph(3=>[1]), false), r2, r4)
# link((a,b)->Link(a,b,BondGraph(2=>[1], 4=>[3]), false), r4, r5)
# link((a,b)->Link(a,b,BondGraph(2=>[1], 4=>[3]), false), r1, r6)
