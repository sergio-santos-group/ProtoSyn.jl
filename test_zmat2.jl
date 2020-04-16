module Treta
# include("src/ProtoSyn.jl")
include("zmat2.jl")

# using .ProtoSyn
using .Zmat

using LinearAlgebra

#lib = ProtoSyn.Peptides.load(ProtoSyn.Peptides.resource_dir)
#mol, state = ProtoSyn.Peptides.build("AGP", lib)


# function res2res(r::ProtoSyn.Residue)
#     res = Residue(r.name, r.id)
#     for at in ProtoSyn.iterbyatom(r)
#         atom = Atom(at.name, at.id, -1, at.symbol)
#         push!(res, atom)
#     end
#     for (pivot, others) in r.bonds
#         append!(res.atoms[pivot].bonds, res.atoms[others])
#     end

#     pivot = get(res, "N")
#     pivot.node.visited = true
#     queue = Atom[pivot]

#     while !isempty(queue)
#         pivot = popfirst!(queue)
#         for atom in pivot.bonds
#             if !atom.node.visited
#                 push!(pivot.node, atom.node)
#                 atom.node.visited = true
#                 push!(queue, atom)
#             end
#         end
#     end
    
#     res
# end

# function mol2top(m::ProtoSyn.Molecule)
#     top = Topology("top", 1, [])#, [])#, [], [])
#     seg = Segment("A", 1)
#     push!(top, seg)
#     res = Residue("", -1)
#     # res = Residue("", -1, seg)
#     id2atom = Dict{Int,Atom}()

#     for (index,at) in enumerate(ProtoSyn.iterbyatom(m))
#         if res.name != at.residue.name || res.id != at.residue.id
#             res = Residue(at.residue.name, at.residue.id)
#             push!(seg, res)
#         end
#         atom = Atom(at.name, at.id, index, at.symbol)
#         id2atom[at.id+at.residue.offset] = atom
#         push!(res, atom)
#         #push!(top, atom)
#     end
#     for (pivot, others) in m.bonds
#         id2atom[pivot].bonds = map(i->id2atom[i], others)
#     end
#     top
# end

# using Base.Cartesian

# function embed(top::Topology, state::State)
#     # assert top.id==state.id
    
#     xyz = zeros(3, state.size)
#     queue = AtomGraphNode[]
#     root = get(ref, "O").node
#     append!(queue, root.children)
    
#     while !isempty(queue)
#         node = popfirst!(queue)
#         append!(queue, node.children)
        
#         i = node.item.index
#         j = node.parent.item.index
#         k = node.parent.parent.item.index
#         istate = state[i]
#         jstate = state[j]        
#         kstate = state[k]        
#         Ri = istate.r
#         Rj = jstate.r
        
#         # local coord system
#         b = istate.b
#         sθ,cθ = sincos(istate.θ)  # angle
#         sϕ,cϕ = sincos(istate.ϕ)  # dihedral
#         x_1 = -b*cθ
#         x_2 =  b*cϕ*sθ
#         x_3 =  b*sϕ*sθ
        
#         # rotate to parent coord system
#         @nexprs 3 u -> vji_u = Rj[u,1]*x_1 + Rj[u,2]*x_2 + Rj[u,3]*x_3
        
#         # UPDATE ROTATION MATRIX
#         @nexprs 3 u -> vjk_u = kstate.t[u] - jstate.t[u]
        
#         # column 1 (x)
#         @nexprs 3 u -> Ri[u,1] = vji_u/b
            
#         # column 3 (z)
#         Zmat.@cross u n_u vji_u vjk_u
#         dn = sqrt(Zmat.@dot(u, n_u, n_u))
#         @nexprs 3 u -> Ri[u,3] = n_u/dn
    
#         # column 2 (y)
#         Zmat.@cross u Ri[u,2] Ri[u,3] Ri[u,1]
        
#         # move to new position
#         @nexprs 3 u -> istate.t[u] = vji_u + jstate.t[u]
#         xyz[:, i] .= istate.t
#     end
#     xyz
# end

# function embed(top::Topology, state::Matrix{Float64})
#     root = get(top.segments[1].residues[1], "N").node
    
#     natoms = length(top)+3
#     xyz = zeros(3, natoms)
#     mat = zeros(3, 3, natoms)
#     @nexprs 3 u -> mat[u,u,root.item.index] = 1
#     queue = AtomGraphNode[]
#     for n in root.children
#         push!(queue, n)
#     end
    
#     while !isempty(queue)
#         node = popfirst!(queue)
#         for n in node.children
#             push!(queue, n)
#         end
#         i = node.item.index

#         r = state[1,i]              # bond
#         sθ,cθ = sincos(state[2,i])  # angle
#         sϕ,cϕ = sincos(state[3,i])  # dihedral
#         x_1 = -r*cθ
#         x_2 = r*cϕ*sθ
#         x_3 = r*sϕ*sθ
        
#         j = node.parent.item.index
        
#         # vector j->i
#         @nexprs 3 u -> vji_u = mat[u,1,j]*x_1 + mat[u,2,j]*x_2 + mat[u,3,j]*x_3
#         println([vji_1 vji_2 vji_3], mat[:,:,j])
        
#         if hasparent(node.parent)
#             k = node.parent.parent.item.index
#             @nexprs 3 u -> vjk_u = xyz[u,k] - xyz[u,j]

#             # column 1 (x)
#             @nexprs 3 u -> mat[u,1,i] = vji_u/r
            
#             # column 3 (z)
#             Zmat.@cross u n_u vji_u vjk_u
#             dn = sqrt(Zmat.@dot(u, n_u, n_u))
#             @nexprs 3 u -> mat[u,3,i] = n_u/dn
        
#             # column 2 (y)
#             Zmat.@cross u mat[u,2,i] mat[u,3,i] mat[u,1,i]
            
#         else
#             mat[:,:,i] .= 0
#             @nexprs 3 u -> mat[u,u,i] = 1
#         end
#         @nexprs 3 u -> xyz[u,i] = vji_u + xyz[u,j]
#     end
#     xyz
# end

# function top2zmat(top::Topology, state::ProtoSyn.State)
    
#     byatom = eachatom(top)
#     natoms = length(byatom)
    
#     zstate = State(natoms)
#     for atom in byatom
#         i = atom.index
#         zstate[i].t .= state.coords[:,i]
#     end
    
#     vij = zeros(3)
#     vjk = zeros(3)
#     vkl = zeros(3)

#     for atom in byatom
#         i = atom.index
#         node = atom.node
        
#         # bond length
#         node = node.parent
#         j = node.item.index
#         @. vij = zstate[j].t - zstate[i].t
#         dij = sqrt(dot(vij,vij))
#         zstate[i].b = dij

#         # angle
#         node = node.parent
#         k = node.item.index
#         @. vjk = zstate[k].t - zstate[j].t
#         djk = sqrt(dot(vjk,vjk))
#         zstate[i].θ = pi - acos(dot(vij,vjk) / (dij*djk))

#         # dihedral
#         node = node.parent
#         l = node.item.index

#         @. vkl = zstate[l].t - zstate[k].t
#         n = cross(vij, vjk)
#         m = cross(vjk, vkl)
#         o = cross(n, m)
#         x = dot(o,vjk)/sqrt(dot(vjk,vjk))
#         y = dot(n,m)
#         zstate[i].ϕ = atan(x,y)




#     end
#     zstate
# end




# function build_tree!(top::Topology)
#     queue = Atom[]
#     iterator = eachatom(top)
#     #println(length(iterator))
#     #println(size(iterator))

#     for atom in iterator
#         atom.node.visited && continue
        
#         atom.node.visited = true
#         push!(queue, atom)
        
#         # i = findfirst(at->!at.node.visited, atom.bonds)
#         # at2 = atom.bonds[i]
#         # at2.node.visited = true
#         # push!(queue, at2)
#         # push!(atom.node, at2)


#         while !isempty(queue)
#             parent = popfirst!(queue)
#             for at in parent.bonds
#                 at.node.visited && continue
#                 push!(parent.node, at.node)
#                 at.node.visited = true
#                 push!(queue, at)
#             end
#         end
#     end
#     root = get(ref, "O").node
#     for (index,atom) in enumerate(iterator)
#         atom.index = index
#         anode = atom.node
#         if !hasparent(anode)
#             push!(root, anode)
#         #    anode.parent = root
#            continue
#         end
#         # !hasparent(anode) && continue
#         child = anode.item.residue
#         parent = anode.parent.item.residue
#         if child===parent || (child.node.visited && parent.node.visited)
#             continue
#         end
#         push!(parent.node, child.node)
#         child.node.visited = true
#         parent.node.visited = true
#     end

#     # head = tail = 0
    
#     # # build atom tree
#     # for atom in top.atoms
        
#     #     atom.node.visited && continue

#     #     #push!(top, atom.node)
#     #     atom.node.visited = true
#     #     tail += 1

#     #     while head < tail
#     #         parent = top.anodes[ head+=1 ]
#     #         for atom in parent.item.bonds
#     #             atom.node.visited && continue
#     #             push!(parent, atom.node)
#     #             #push!(top, atom.node)
#     #             atom.node.visited = true
#     #             tail += 1
#     #         end
#     #     end

#     # end

#     # build residue tree
#     #for residue in iterbyresidue(top)
#     #    push!(top.rnodes, ResidueGraphNode(residue, nothing))
#     #end
    
#     # for anode in top.anodes
#     #     (anode.parent === nothing) && continue
#     #     child = anode.item.residue
#     #     parent = anode.parent.item.residue
#     #     if child===parent || (child.node.visited && parent.node.visited)
#     #         continue
#     #     end
#     #     # child.node.parent = parent.node
#     #     # push!(parent.node.children, child.node)
#     #     push!(parent.node, child.node)
#     #     child.node.visited = true
#     #     parent.node.visited = true
#     # end

# end

top, state = read("resources/Peptides/old/ala.pdb", Zmat.PDB)
build_tree!(top)
sync!(state, top)

#top = mol2top(mol)
#build_tree!(top)
#zmat = Treta.top2zmat(Treta.top, Treta.state)
#xyz = Treta.embed(Treta.top,zmat)

# Base.merge!(combine::Function, at1::Atom, at2::Atom) = begin
#     top1 = at1.residue.segment.topology
#     top2 = at2.residue.segment.topology
#     if top1 === top2
#     else
#     end
# end

join(r1::Residue, r2::Residue) = begin
    
    atC = get(r1, "C")
    atN = get(r1, "N")
    
    push!(atC.bonds, atN)
    push!(atN.bonds, atC)

    push!(r1.node, r2.node)
    push!(atC.node, atN.node)

end

# peptide specific
insert!(r1::Opt{Residue}, r2::Opt{Residue}, rn::Residue) = begin
    if r1===r2
        error("r1 and r2 must be different")
    elseif rn===r1 || rn===r2
        error("rn cannot be equal to r1 or r2")
    elseif !isorphan(rn)
        error("the residue to be inserted must be orphan")
    end
    
    if r1!==nothing && r2 !==nothing
        split(r1, r2)
        join(r1, rn)
        join(rn, r2)
    elseif r1===nothing
        join(rn, r2)
    else
        join(r1, rn)
    end
end

split(r1::Residue, r2::Residue) = begin
    # split the residue graph
    split(r1.node, r2.node)

    # break the C-N amide bond
    atC = get(r1, "C")
    atN = get(r2, "N")
    i = findfirst(at->at===atN, atC.bonds)
    i !== nothing && deleteat!(atC.bonds, i)
    
    i = findfirst(at->at===atC, atN.bonds)
    i !== nothing && deleteat!(atN.bonds, i)

    # and split the atom graph
    split(atC.node, atN.node)

end


# using Printf

# Base.write(io::IO, top::Topology, xyz::Matrix{Float64}, ::Type{Zmat.PDB}) = begin
    
#     println(io, "MODEL")
#     node = get(top.segments[1].residues[1], "N").node

#     # traverse(node) do anode
#     #     atom = anode.item
#     for atom in eachatom(top)
#         s = @sprintf("ATOM %6d %4s %-4sA %3d    %8.3f%8.3f%8.3f%24s",
#             atom.index, atom.name,
#             atom.residue.name, atom.residue.id,
#             xyz[1,atom.index],
#             xyz[2,atom.index],
#             xyz[3,atom.index],
#             atom.symbol)
#         println(io, s)
#     end

#     # traverse(node) do anode
#     for atom in eachatom(top)
#         print(io, @sprintf("CONECT%5d", atom.index))
#         foreach(n->print(io, @sprintf("%5d",n.item.index)), atom.node.children)
#         println(io,"")
#     end
#     println(io, "ENDMDL")
# end

# Base.write(io::IO, top::Topology, state::ProtoSyn.State, ::Type{ProtoSyn.PDB}) = begin
#     Base.write(io, top, state.coords, ProtoSyn.PDB)
# end


end
