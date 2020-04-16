using ..ProtoSyn
using LinearAlgebra: det, svd!

@doc """
brute force peptide link.
A new link C->N is always created if atoms named "C" and "N" exist
in residues lr1 and lr2, respectively.
"""
function peptidebond(r1::Residue, r2::Residue)
    if (atC = get(r1, "C", nothing)) === nothing ||
       (atN = get(r2, "N", nothing)) === nothing
       return nothing
    end
    Link(r1, r2, atC.id => atN.id)
end

function fit_using_bb(
    old::Residue, xold::AbstractArray,
    new::Residue, xnew::AbstractArray, state::State)

    iold = [
        get(old,  "N").id
        get(old, "CA").id
        get(old,  "C").id
    ]

    inew = [
        get(new,  "N").id
        get(new, "CA").id
        get(new,  "C").id
    ]
    target = xold[:,iold]
    mobile = xnew[:,inew]

    target_centroid = sum(target, dims=2)/3
    mobile_centroid = sum(mobile, dims=2)/3
    target .-= target_centroid
    mobile .-= mobile_centroid
    
    f = svd!(mobile * target')
    if det(f.Vt) * det(f.U) < 0.0
        f.Vt[3,:] .*= -1.0
    end
    # R = f.Vt' * f.U'
    R = (f.U * f.Vt)'
    # println(f.Vt' * f.U' - (f.U * f.Vt)')
    xnew .-= mobile_centroid
    xnew .= R*xnew .+ target_centroid

    i = get(new, "O")
    j = get(old, "O")
    if (i!==nothing) && (j!==nothing)
        xnew[:,i.id] .= xold[:,j.id]
    end

    if (i = get(new, "H")) !== nothing
        if (j = get(old, "H")) !== nothing
            xnew[:,i.id] .= xold[:,j.id]
        else
            ProtoSyn.cproduct(new, ("-C","N","CA")) do atoms
                C,N,CA = atoms
                xC = state.coords[:,C.id+C.residue.offset]
                xN = xnew[:,N.id]
                xCA = xnew[:,CA.id]
                vNC = xC - xN
                vNCA = xCA - xN
                dNC = norm(vNC)
                dNCA = norm(vNCA)
                P = (dNC*vNCA + dNCA*vNC)/(dNC*dNCA)
                xnew[:,i.id] .= xN - P/norm(P)
            end
        end
    end
end


#function Base.replace!(mol::Molecule, state::State, old_new::Pair{Residue,Residue})
#    replace!(fit_using_bb, mol, state, old_new)
#end
function replace!(mol::Molecule, state::State, old_new::Pair{Residue,Residue})
   Base.replace!(fit_using_bb, mol, state, old_new, peptidebond)
end


function build(::Type{T}, sequence::Vector{String}, lib::ResidueLib) where {T<:AbstractFloat}
    
    mol = Molecule()
    mol.name = "peptide"

    foreach(sequence) do aa
        template = lib[aa]
        @assert ProtoSyn.istemplate(template) "peptides can only be built from template residues"
        push!(mol, deepcopy(template))
    end
    
    seqlen = length(sequence)
    
    # identify root residue. Since a peptide is a linear chain,
    # the root is always identified as the first residue
    if seqlen > 0
        mol.root = mol.residues[1]
    end
    
    # link all aminoacids into a single chain
    for n=1:seqlen-1
        ProtoSyn.link(peptidebond, mol.residues[n], mol.residues[n+1])
    end
    
    for (n,res) in enumerate(mol.residues)
        res.id=n
    end
    ProtoSyn.update!(mol)
    
    # instantiate a state for the recently created molecule
    # and set up coordinates accordingly. Some of the magic numbers
    # below are required to ensure that the C-N bond length is
    # ca. 1.34 Å and the CA-C-N(i+1) and C(n-1)-N-CA angles
    # around 120 degrees. Also, and becuase each residue is no longer a
    # template, the :coords entry is removed from the metadata dictionary.
    state = State{T}(size(mol))
    # state.coords = vcat([
    state.coords = hcat([
       (
           dx=3.67*(r.id-1);
           f=(-1)^(r.id);
        #    x=pop!(r.metadata,:coords);
           (r.coords .+ [dx, -0.28, 0]) .* [1, f, f]
        #    (x .+ [dx -0.28 0]) .* [1 f f]
        )
        for r in mol.residues]...
    )

    return mol,state
end

build(seq::Vector{String}, lib::ResidueLib) = build(Float64, seq, lib)
build(seq::String, lib::ResidueLib) = build([one_2_three[c] for c in seq], lib)

function finddihedrals(mol::Molecule; backbone::Bool=true, sidechain::Bool=true)
    if !isvalid(mol)
        error("finddihedrals requires an updated molecule")
    end

    # if this molecule is not coherent (not yet updated
    # after creation of modification) or no bonded table
    # exists, then return nothing!
    if !isvalid(mol) || isempty(mol.bonds)
        return nothing
    end

    # container for all dihedrals
    dihedrals = Vector{Dihedral}()
    
    # callback function for ProtoSyn.cproduct:
    #  it takes a 4-tuple of indices, instantiates a new
    #  dihedral ofthe given type and add it to the container
    cb(atoms::Vector{Atom}, t::DihedralTypes.Type) = push!(dihedrals, Dihedral(atoms..., t))
    
    # temporary array for storing atoms
    ar = Vector{Atom}(undef, 4)
    
    for res in mol.residues
        if backbone
            ProtoSyn.cproduct(j->cb(j, DihedralTypes.phi),   res, ("-C", "N","CA",  "C"), 0, ar)
            ProtoSyn.cproduct(j->cb(j, DihedralTypes.psi),   res, ( "N","CA", "C", "+N"), 0, ar)
            ProtoSyn.cproduct(j->cb(j, DihedralTypes.omega), res, ("CA", "C","+N","+CA"), 0, ar)
        end
        if sidechain
            ProtoSyn.cproduct(j->cb(j, DihedralTypes.chi1), res, ( "N", "CA", "CB", "CG"), 0, ar)
            ProtoSyn.cproduct(j->cb(j, DihedralTypes.chi2), res, ("CA", "CB", "CG", "CD"), 0, ar)
            ProtoSyn.cproduct(j->cb(j, DihedralTypes.chi3), res, ("CB", "CG", "CD", "CE"), 0, ar)
            ProtoSyn.cproduct(j->cb(j, DihedralTypes.chi4), res, ("CG", "CD", "CE", "CZ"), 0, ar)
            ProtoSyn.cproduct(j->cb(j, DihedralTypes.chi5), res, ("CD", "CE", "CZ", "CH"), 0, ar)
        end
    end

    dihedrals

end


# function findcrankshafts(mol::Molecule)
#     if !isvalid(mol)
#         error("findcrankshafts requires an updated molecule")
#     end
    
# end






function traverse(u::Int, p::Int,
    color::Vector{Int},
    mark::Vector{Int},
    par::Vector{Int},
    cyclenumber::Int,
    bonds::BondGraph)
    if color[u] == 2
        return cyclenumber
    end
    
    if color[u] == 1
        cur = p
        #cyclenumber += 1
        cyclenumber = u
        mark[cur] = cyclenumber
        while cur != u
            cur = par[cur]
            mark[cur] = cyclenumber
        end
        return cyclenumber
    end

    par[u] = p
    color[u] = 1

    for v in bonds[u]
        v == par[u] && continue
        cyclenumber = traverse(v, u, color, mark, par, cyclenumber, bonds)
    end
    color[u] = 2
    return cyclenumber
end

function treta(res::Residue)
    #ts = ProtoSyn.TraverseState(size(res))
    #reset(ts)

    mark  = zeros(Int, size(res))
    color  = zeros(Int, size(res))
    parent = zeros(Int, size(res))
    
    traverse(1, 1, color, mark, parent, 0, res.bonds)
    for i = 1:size(res)
        if length(res.bonds[i]) == 1
            j = res.bonds[i][1]
            if mark[j]>0
                mark[i] = mark[j]
            end
        end

    end

    println(mark)
    println(color)
    println(parent)

    m = Dict{Int, Vector{Int}}()
    # visited = falses(size(res))
    # for i=1:size(res)
    #     length(res.bonds[i]) == 1 && continue
    #     for j in res.bonds[i]
    #         visited[j] && continue
            
    #         #visited[j] = true
    #         push!(get!(m, mark[i] == 0 ? i : mark[i], []), j)

    #     end
    #     visited[i] = true
    # end
    ts = ProtoSyn.TraverseState(size(res))
    
    N = get(res, "N").id
    ts.visited[N] = true
    head = tail = 0
    ts.indices[tail+=1] = N
    while head < tail
        id1 = ts.indices[head+=1]
        
        length(res.bonds[id1]) == 1 && continue

        for id2 in res.bonds[id1]
            if !ts.visited[id2]
                ts.visited[id2] = true
                ts.indices[tail+=1] = id2

                push!(get!(m, mark[id1] == 0 ? id1 : mark[id1], []), id2)

            end
        end
    end
    println(m)
end