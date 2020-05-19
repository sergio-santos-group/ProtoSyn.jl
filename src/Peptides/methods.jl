using LinearAlgebra: norm, cross
using ..ProtoSyn

@doc """
reset aminoacids orientation such that the N is at the origin,
the N-C vector lays along the x-axis and the N-CA-C atoms are
on the xy plane.
"""
function reset!(aminoacids::Dict{String, Residue})
    
    rmat = zeros(3,3)
    
    # for aa in values(aminoacids)
    for name in keys(three_2_one)
        
        aa = get(aminoacids, name, nothing)
        (aa===nothing) && continue
        
        # gather coordinates into a single matrix
        xyz = Array(hcat(([at.x,at.y,at.z] for at in aa.atoms)...)')

        # identify backbone atoms. The residue will be aligned
        # such that the N and C atoms lay along the x-axis and the
        # N, CA, and C atoms in the xy plane
        idN  = findfirst(a->a.name=="N",  aa.atoms)
        idCA = findfirst(a->a.name=="CA", aa.atoms)
        idC  = findfirst(a->a.name=="C",  aa.atoms)
        
        # shift coordinates so that atom N is at origin
        xyz .-= xyz[idN, :]'
        
        # align the N->C vector along x-axis
        theta = acos(xyz[idC,1] / norm(xyz[idC,:]))
        axis = cross([1.0, 0.0, 0.0], xyz[idC,:])
        rotmat!(rmat, axis, -theta)
        xyz = xyz * rmat'
        
        # put the N->CA vector on the xy-plane
        v = xyz[idCA,:]
        theta = sign(v[3])*acos(v[2] / sqrt(v[2]^2 + v[3]^2))
        rotmat!(rmat, [1.0, 0.0, 0.0], -theta)
        xyz = xyz * rmat'

        # update residue coordinates
        for i =1:length(aa.atoms)
            aa.atoms[i].x = xyz[i,1]
            aa.atoms[i].y = xyz[i,2]
            aa.atoms[i].z = xyz[i,3]
        end
    end

    aminoacids
end

# function peptidelink(lr1::LinkedResidue, lr2::LinkedResidue)
    
#     # the peptide bond is C[n]->N[n+1]. If either the
#     # C or the N is not found, then a peptide bond
#     # cannot be made.

#     # the peptide bond can be made either via
#     # C(lr1)->N(lr2) or C(lr2)->N(lr1) because
#     # the order by which both residues are given
#     # might not reflect their ordering within the
#     # chain. Hence, both alternatives are possible.
#     # Nevertheless, the input ordering is first
#     # tried and, only if not possible to form a bond,
#     # will the second alternative be tried.
#     atC = get(lr1, "C", nothing)
#     atC == nothing && return
#     idC = atC.id
    
#     atN = get(lr2, "N", nothing)
#     atN == nothing && return
#     idN = atN.id

#     # both indices have been found. We now check
#     # if both accept an additional bond. These are the
#     # intra residue bonds each atom has:
#     degC = length(lr1.source.bonds[idC])
#     degN = length(lr2.source.bonds[idN])

#     # now we update for possible additional
#     # inter residue bonds
#     degC += count(l->l.residue1 === lr1 && haskey(idC, l.bonds), lr1.links)
#     degC > 2 && return

#     degN += count(l->l.residue2 === lr1 && haskey(idN, l.bonds), lr2.links)
#     degN > 2 && return

#     Link(lr1, lr2, ConnectGraph(idC => [idN]))
# end


@doc """
brute force peptide link.
A new link C->N is always created if atoms named "C" and "N" exist
in residues lr1 and lr2, respectively.
"""
function peptidebond(lr1::LinkedResidue, lr2::LinkedResidue)
    if (atC = get(lr1, "C", nothing)) === nothing ||
       (atN = get(lr2, "N", nothing)) === nothing
       return nothing
    end
    Link(lr1, lr2, ConnectGraph(atC.id => [atN.id]))
end



function build(::Type{T}, letters::String, aminoacids::ResidueLib) where {T}
    build(T, [one_2_three[letter] for letter in letters], aminoacids)
end

build(letters::String, aminoacids::ResidueLib) = build(Float64, letters, aminoacids)

build(sequence::Vector{String}, aminoacids::ResidueLib) = build(Float64, sequence, aminoacids)

function build(::Type{T}, sequence::Vector{String}, aminoacids::ResidueLib) where {T}
    
    # reset aminoacid orientation
    reset!(aminoacids)

    # instantiate a new molecule bearing the requested aminoacids
    mol = Molecule(
        name="peptide",
        residues = [
            LinkedResidue(id=id, source=aminoacids[name])
            for (id,name) in enumerate(sequence)
        ]
    )
    

    # generate links for all consecutively bonded residues
    for n = 1:length(sequence)-1
        # get a pair of consecutive l-residues (aminoacids)
        # and attempt at linking them. if successfull,
        # add that link to the molecule.
        link = ProtoSyn.bind(peptidebond, mol.residues[n], mol.residues[n+1], mol)
        
        # push this link (if viable) to the molecule (maybe this should
        # be done by the core bind function?)
        # link!==nothing && push!(mol, link)
        
    end

    
    # update molecule (this steps also determines the
    # size of the molecule)
    ProtoSyn.update!(mol)

    # instantiate a state for the recently created molecule
    # and set up coordinates accordingly. Some of the magic numbers
    # below are required to ensure that the C-N bond length is
    # ca. 1.34 Å and the CA-C-N(i+1) and C(n-1)-N-CA angles
    # around 120 degrees.
    state = State{T}(mol.size)
    for (n, lr) in enumerate(mol.residues)
        # the amount to move along the x axis
        δx = 3.670*(n-1)
        # δx = 0.3670*(n-1)
        # rotation factor for odd residues (rotation
        # of π about the x-axis)
        f = (-1)^n
        for i=1:length(lr)
            atom = lr.source.atoms[i]
            i += lr.offset
            state.coords[i,1] = atom.x + δx
            state.coords[i,2] = f*(atom.y - 0.28)
            # state.coords[i,2] = f*(atom.y - 0.028)
            state.coords[i,3] = f*atom.z
        end
    end

    # return molecule & state
    mol, state
end


function setss!(state::State, mol::Molecule, rng::UnitRange{Int}, conf::Symbol)
    
    if conf == :antiparallel_sheet
        ϕ, ψ = deg2rad(-139.0), deg2rad(135.0)
    elseif conf == :parallel_sheet
        ϕ, ψ = deg2rad(-119.0), deg2rad(113.0)
    elseif conf == :helix
        ϕ, ψ = deg2rad(-60.0), deg2rad(-45.0)
    else
        return state
    end

    ProtoSyn.set!(state, mol, rng, ("-C","N","CA","C"), ϕ)
    ProtoSyn.set!(state, mol, rng, ("N","CA","C","+N"), ψ)
    state
end



# function finddihedrals(mol::Molecule;
#     backbone::Bool=true, sidechain::Bool=true)

#     # if this molecule is not coherent (not yet updated
#     # after creation of modification) or no bonded table
#     # exists, then return nothing! 
#     if !mol.coherent || (mol.bonds === nothing)
#         return nothing
#     end

#     # container for all dihedrals
#     dihedrals = Vector{Dihedral}()

#     # names of all atoms
#     atnames = map(at->at.name, ProtoSyn.iterbyatom(mol))

#     # auxilliary function that finds the first occurrence
#     # of an atom named "name2" in the list of all atoms
#     # bonded to "pivot". If found, it assigns "dtype" to the
#     # generated dihedral.
#     function findpair(i2::Int, s1::String, s3::String, s4::String, dtype::DihedralTypes.Type)
#         # i1--i2--i3--i4
#         i1 = findfirst(n-> atnames[n] == s1, mol.bonds[i2])
#         i1 === nothing && return
        
#         i3 = findfirst(n-> atnames[n]==s3, mol.bonds[i2])
#         i3 === nothing && return
#         i3 = mol.bonds[i2][i3]
        
#         i4 = findfirst(n-> atnames[n]==s4, mol.bonds[i3])
#         i4 === nothing && return 
        
#         push!(dihedrals,
#             Dihedral(
#                 a0=mol.bonds[i2][i1],
#                 a1=i2,
#                 a2=i3,
#                 a3=mol.bonds[i3][i4],
#                 start=i3,
#                 type=dtype)
#         )
#     end


#     for (id1, atom) in enumerate(ProtoSyn.iterbyatom(mol))    
        
#         if backbone
#             if atom.name=="N"
#                 #  ϕ  (phi)
#                 findpair(id1, "C", "CA", "C", DihedralTypes.phi)
#             elseif atom.name=="C"
#                 #  ω  (omega)
#                 findpair(id1, "CA", "N", "CA",  DihedralTypes.omega)
#             elseif atom.name=="CA"
#                 #  ψ  (psi)
#                 findpair(id1, "N", "C", "N",  DihedralTypes.psi)
#             end
#         end

#         if sidechain
#             if atom.name=="CA"
#                 # χ1  (chi1)
#                 findpair(id1, "N", "CB", "CG", DihedralTypes.chi1)
#             elseif atom.name=="CB"
#                 # χ2  (chi2)
#                 findpair(id1, "CA", "CG", "CD", DihedralTypes.chi2)
#             elseif atom.name=="CG"
#                 # χ3  (chi3)
#                 findpair(id1, "CB", "CD", "CE", DihedralTypes.chi3)
#             elseif atom.name=="CD"
#                 # χ4  (chi4)
#                 findpair(id1, "CG" ,"CE", "CZ", DihedralTypes.chi4)
#             elseif atom.name=="CE"
#                 # χ5  (chi5)
#                 findpair(id1, "CD", "CZ", "CH", DihedralTypes.chi5)
#             end
#         end

#     end

#     dihedrals

# end

function finddihedrals(mol::Molecule; backbone::Bool=true, sidechain::Bool=true)

    # if this molecule is not coherent (not yet updated
    # after creation of modification) or no bonded table
    # exists, then return nothing!
    if !isvalid(mol) || (mol.bonds === nothing)
        return nothing
    end

    # container for all dihedrals
    dihedrals = Vector{Dihedral}()
    
    # callback function for ProtoSyn.cproduct:
    #  it takes a 4-tuple of indices, instantiates a new
    #  dihedral ofthe given type and add it to the container
    apply(i::Vector{Int}, t::DihedralTypes.Type) = push!(dihedrals, Dihedral(i..., t))
    
    # temporary array for storing atom indices
    ar = zeros(Int, 4)
    
    for res in mol.residues
        if backbone
            ProtoSyn.cproduct(j->apply(j, DihedralTypes.phi),   res, ("-C", "N","CA",  "C"), 0, ar)
            ProtoSyn.cproduct(j->apply(j, DihedralTypes.psi),   res, ( "N","CA", "C", "+N"), 0, ar)
            ProtoSyn.cproduct(j->apply(j, DihedralTypes.omega), res, ("CA", "C","+N","+CA"), 0, ar)
        end
        if sidechain
            ProtoSyn.cproduct(j->apply(j, DihedralTypes.chi1), res, ( "N", "CA", "CB", "CG"), 0, ar)
            ProtoSyn.cproduct(j->apply(j, DihedralTypes.chi2), res, ("CA", "CB", "CG", "CD"), 0, ar)
            ProtoSyn.cproduct(j->apply(j, DihedralTypes.chi3), res, ("CB", "CG", "CD", "CE"), 0, ar)
            ProtoSyn.cproduct(j->apply(j, DihedralTypes.chi4), res, ("CG", "CD", "CE", "CZ"), 0, ar)
            ProtoSyn.cproduct(j->apply(j, DihedralTypes.chi5), res, ("CD", "CE", "CZ", "CH"), 0, ar)
        end
    end

    dihedrals

end


function findcrankshafts(mol)
    
    # pointer to connectivity graph
    graph = mol.bonds
    
    visited = falses(mol.size)
    
    # id-to-name map, restricted to backbone atoms
    atnames = Dict(
        i => atom.name
        for (i,atom) in enumerate(ProtoSyn.iterbyatom(mol))
        if atom.name ∈ ("N", "CA", "C")
    )
    
    # mark all atoms in flagged residues as true.
    # This implies that all flagged atoms wil not be
    # able to define a crankshaft.
    flagged = falses(mol.size)
    for lr in mol.residues
        if lr.flag
            for atom in lr.source.atoms
                flagged[lr.offset+atom.id] = true
            end
        end
    end

    crankshafts = AxisRotatableBlock[]

    function traverse(pivot::Int, cur::Int, prev::Int)
        # is this a CA atom and not the pivot? if so
        # this constitutes as possible crankshaft pair and
        # one only has to identify the ordering because we
        # assume that travel direction is restricted
        # to CA->C->N->(...)->N->CA.
        # Because peptides are linear chains, and the full
        # traversal is restricted to the backbone atoms, then
        # the last CA atom has to be reached through "N".
        # If not, we simply reverse pair ordering.
        if !flagged[cur] && (cur > pivot) && (atnames[cur] == "CA")
            if atnames[prev] != "N"
                cur, pivot = pivot, cur
            end
            
            # identify the ID of the first atom for future traversal
            idC = findfirst(i->get(atnames, i, nothing) == "C", graph[pivot])
            if idC !== nothing
                push!(
                    crankshafts,
                    AxisRotatableBlock(pivot,cur,graph[pivot][idC])
                )
            end
        end

        # mark the current index as visited
        # and do the traversal
        visited[cur] = true
        for i in graph[cur]
            !visited[i] && haskey(atnames,i) && traverse(pivot, i, cur)
        end
        visited[cur] = false
    end


    # find all candidates
    for i in keys(atnames)
        if atnames[i] == "CA" && !flagged[i]
            traverse(i, i, i)
        end
    end

    return crankshafts

end
