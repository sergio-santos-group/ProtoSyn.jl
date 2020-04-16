using Printf: @sprintf

function repr_opt(item::AbstractArray)
    "$(size(item,1))-element array"
end


function repr_opt(item::AbstractDict)
    "$(length(item))-element dict"
end


function repr_opt(item::Nothing)
    "n.d."
end


#region Base.get
# -----------------------------------------------

Base.get(lr::LinkedResidue, k::AbstractString, default) = begin
    get(lr.source.atomsbyname, k, default)
end
Base.get(lr::LinkedResidue, k::AbstractString) = get(lr, k, nothing)



Base.get(r::Residue, k::AbstractString, default) = begin
    get(r.atomsbyname, k, default)
end
Base.get(r::Residue, k::AbstractString, default) = get(r, k, nothing)

#endregion Base.get


#region Base.length
# -----------------------------------------------

Base.length(r::Residue) = length(r.atoms)

Base.length(lr::LinkedResidue) = length(lr.source)

#endregion Base.length



#region Base.show
# -----------------------------------------------

Base.show(io::IO, item::Atom) = begin
    print(io, "Atom:")
    print(io, "\n     id = $(item.id)")
    print(io, "\n   name = $(item.name)")
    print(io, "\n symbol = $(item.symbol)")
    print(io, "\n parent = $(item.parent === nothing ? "n.d." : item.parent.name)")
end


Base.show(io::IO, item::Residue) = begin
    print(io, "$(length(item))-atom Residue:")
    print(io, "\n   name = $(item.name)")
    print(io, "\n  atoms = $(repr_opt(item.atoms))")
    print(io, "\n  bonds = $(repr_opt(item.bonds))")
end


Base.show(io::IO, lr::LinkedResidue) = begin
    print(io, "$(length(lr))-atom LinkedResidue:")
    print(io, "\n      id = $(lr.id)")
    print(io, "\n  source = $(lr.source.name)")
    print(io, "\n  offset = $(lr.offset)")
    print(io, "\n   links = $(repr_opt(lr.links))")
end


Base.show(io::IO, link::Link) = begin
    print(io, "$(length(link.bonds))-bond Link:")
    print(io, "\n  residue1 = $(link.residue1.source.name).$(link.residue1.id)")
    print(io, "\n  residue2 = $(link.residue2.source.name).$(link.residue2.id)")
    if !isa(link.bonds, Nothing)
        print(io, "\n  bonds:")
        for pair in pairs(link.bonds)
            print(io, "\n   $(pair.first) => $(join(map(string, pair.second), ", "))")
        end
    end
end


Base.show(io::IO, mol::Molecule) = begin
    for (rid,lr) in enumerate(mol.residues)
        offset = mol.offset + lr.offset
        for at in lr.source.atoms
            atid = at.id + offset
            s = @sprintf("ATOM %6d %4s %-4s %3d %2s",
                atid, at.name,
                lr.source.name, rid, at.symbol)
            s2 = map(x->@sprintf("%5d", x+mol.offset), mol.bonds[atid])
            println(io, s, " -> ", mol.bonds[atid] .+ mol.offset)
        end
    end
end

#endregion Base.show



#region Base.push!
# -----------------------------------------------

Base.push!(residue::Residue, atom::Atom) = begin
    
    #if residue.atoms === nothing
    #    residue.atoms = Atom[]
    #    residue.atomsbyname = Dict{String, Atom}()
    #end

    if haskey(residue.atomsbyname, atom.name)
        error("Atom names within residues must be unique!")
    end
    
    if atom.parent !== nothing
        error("Atom already has a parent!")
    end

    atom.parent = residue
    atom.id = length(residue) + 1
    
    push!(residue.atoms, atom)
    residue.atomsbyname[atom.name] = atom

    residue
end


Base.push!(mol::Molecule, lr::LinkedResidue) = begin
    lr.offset = mapreduce(r->length(r), +, mol.residues; init=0)
    lr.id = size(mol.residues, 1) + 1
    push!(mol.residues, lr)
    invalidate!(mol)
    mol
end


Base.push!(mol::Molecule, link::Link) = begin
    push!(mol.links, link)
    invalidate!(mol)
    mol
end


Base.push!(lresidue::LinkedResidue, link::Link) = begin
    push!(lresidue.links, link)
    lresidue
end
#endregion


function invalidate!(mol::Molecule)
    mol.coherent = false
    mol
end



Base.isvalid(mol::Molecule) = mol.coherent

# function Base.delete!(mol::Molecule, fragment::LinkedResidue)
#     # find location of the given fragment in the molecule
#     idx = findfirst(f->f===fragment, mol.fragments)
    
#     # if the given linked fragment is not found in this
#     # molecule, then do nothing
#     if idx === nothing
#         return mol
#     end

#     # remove the fragment
#     deleteat!(mol.fragments, idx)
    
#     # remove all links this fragment might have from
#     # from the molecule's link list
#     for link in fragment.links
#         idx = findfirst(l->l===link, mol.links)
#         if idx !== nothing
#             deleteat!(mol.links, idx)
#         end
#         # remove this link from both fragments
#         if link.residue1 !== fragment
#             delete!(link.residue1, link)
#         else
#             delete!(link.residue2, link)
#         end
#     end
#     # mark this molecule as no longer being coherent 
#     mol.coherent = false
#     mol
# end


# function Base.delete!(frag::LinkedResidue, link::Link)
#     idx = findfirst(l->l===link, frag.links)
#     if idx !== nothing
#         deleteat!(frag.links, idx)
#     end
#     frag
# end


@doc """


finds atom indices corresponding to the requested atom names, based on
the given residue. Names starting with '-'/'+' are assumed to belong to the
previous/next residue. If multiple residues are connected to the central (given)
residue, all tuple are generated.
"""
function cproduct(f::Function, residue::LinkedResidue, names::Tuple{Vararg{String}})
    _cproduct(f, residue, names, 0, zeros(Int, length(names)), residue)
end

function cproduct(f::Function,
    residue::LinkedResidue, names::Tuple{Vararg{String}}, n::Int, state::Vector{Int})
    _cproduct(f, residue,names,n,state,residue)
end

function _cproduct(f::Function,
    residue::LinkedResidue, names::Tuple{Vararg{String}}, n::Int, state::Vector{Int}, pres::LinkedResidue)
    
    # apply function if the current state is complete
    if n==length(names)
        f(state)
    else
        name = names[n+1]
        if startswith(name, '+')
            for link in residue.links
                # is link.residue1 the source (central) residue? if so, the
                # next ('+') residue is residue2.
                # It is assumed that every link is directed from
                # residue1 (source) -> residue2 (target)
                if link.residue1 === residue
                    (n>0) && startswith(names[n],'+') && (link.residue2 !== pres) && continue
                    atom = get(link.residue2, SubString(name,2), nothing)
                    (atom===nothing) && continue
                    state[n+1] = atom.id + link.residue2.offset
                    _cproduct(f, residue, names, n+1, state, link.residue2)
                end
            end
        elseif startswith(name, '-')
            for link in residue.links
                # is link.residue2 the target residue? if so, the
                # previous ('-') residue is residue1
                if link.residue2 === residue
                    (n>0) && startswith(names[n],'-') && (link.residue1 !== pres) && continue
                    atom = get(link.residue1, SubString(name,2), nothing)
                    (atom===nothing) && continue
                    state[n+1] = atom.id + link.residue1.offset
                    _cproduct(f, residue, names, n+1, state, link.residue1)
                end
            end
        else
            atom = get(residue, name, nothing)
            if atom!==nothing
                state[n+1] = atom.id + residue.offset
                _cproduct(f, residue, names, n+1, state, residue)
            end
        end
    end

end







export set!

function set!(state::State, mol::Molecule, rng::UnitRange{Int}, atnames::NTuple{4, String}, θ::Float64)
    
    xyz = state.coords

    # Allocate-once variables
    rmat = zeros(3,3)           # rotation matrix
    mask = falses(mol.size)     # mask for graph traversal
    indices = zeros(Int, 4)     # container

    for i in rng
        cproduct(mol.residues[i], atnames, 0, indices) do idxs
            dihd = Dihedral(idxs...)
            rotate!(xyz, mol.bonds, dihd, θ-measure(dihd, xyz), mask, rmat)
        end
    end
    state
end

function set!(state::State, mol::Molecule, dihedral::Dihedral, θ::Float64)
    rotate!(state, mol, dihedral, θ-measure(dihedral, state.coords))
end

# function set!(state::State, mol::Molecule, rng::UnitRange{Int}, atnames::NTuple{3, String}, θ::Float64)
    
#     xyz = state.coords

#     # Allocate-once variables
#     rmat = zeros(3,3)           # rotation matrix
#     mask = falses(mol.size)     # mask for graph traversal
#     indices = zeros(Int, 3)     # container
    
#     foreach(mol.residues) do residue
#         cproduct(residue, atnames, 0, indices) do idxs
#             dihd = Angle(idxs...)
#             rotate!(xyz, mol.bonds, dihd, θ-measure(dihd, xyz), mask, rmat)
#         end
#     end


# end