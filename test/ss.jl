module Tst

include("../src/ProtoSyn.jl")

using .ProtoSyn
using .ProtoSyn.Peptides


# construct model
reslib = include(joinpath(ProtoSyn.Peptides.resource_dir, "aminoacids.jl"))
mol, state = ProtoSyn.Peptides.build("?"^10, reslib)



# function work(N)
#     for i=1:N
#         ProtoSyn.Peptides.setss!(state, mol, 10:30, :helix)
#         ProtoSyn.Peptides.setss!(state, mol, 40:60, :helix)
#         ProtoSyn.Peptides.setss!(state, mol, 70:95, :sheet)
#     end
# end

# # @time work(1)
# # @time work(1e2)



# dihedrals = ProtoSyn.Peptides.finddihedrals(mol; sidechain=false)

# phis = filter(d->d.type==ProtoSyn.Peptides.DihedralTypes.phi, dihedrals)
# psis = filter(d->d.type==ProtoSyn.Peptides.DihedralTypes.psi, dihedrals)


block1 = 2:4
block2 = 6:9
#block3 = 70:95

# for block in (block1,block2)#block2,block3)
#     for i in block
#         mol.residues[i].flag = true
#     end
# end

#foreach(lr->lr.flag=false, mol.residues)

crankshafts = ProtoSyn.Peptides.findcrankshafts(mol)




fout = open("../tmp.pdb", "w")
write(fout, mol, state)
# ProtoSyn.Peptides.setss!(state, mol, block1, :helix)
# write(fout, mol, state)
# ProtoSyn.Peptides.setss!(state, mol, block2, :helix)
# write(fout, mol, state)
# ProtoSyn.Peptides.setss!(state, mol, block3, :sheet)
# write(fout, mol, state)


# crankshafts2 = filter(crankshafts) do cr
#     for block in (block1,block2,block3)

#         r1 = mol.residues[block.start]
#         at1 = get(r1, "CA", nothing)
#         id1 = at1.id + r1.offset

#         r2 = mol.residues[block.start]
#         at2 = get(r2, "CA", nothing)
#         id2 = at2.id + r2.offset
        
#         if (id1 <= cr.a1 <= id2) || (id1 <= cr.a2 <= id2)
#             return false
#         end
#     end
#     true
# end


# for phi in phis
#     ProtoSyn.set!(state, mol, phi, deg2rad(90.0))
# end

write(fout, mol, state)



close(fout)



# dihedrals = ProtoSyn.Peptides.finddihedrals(mol; sidechain=false)

end