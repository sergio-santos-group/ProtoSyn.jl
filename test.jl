module Tst


include("src/ProtoSyn.jl")
using .ProtoSyn
using .ProtoSyn.Forcefields

const PS = ProtoSyn



# # frags = PS.loadresidues(joinpath(PS.resource_dir, "Peptides", "aminoacids.pdb"))
# residues = ProtoSyn.Peptides.loadresidues()

# # ff = ProtoSyn.Forcefields.load("amber03")
ff = include(joinpath(PS.Forcefields.resource_dir, "amber03/forcefield.jl"))

aa = include(joinpath(PS.Peptides.resource_dir, "aminoacids.jl"))

mol,state = ProtoSyn.Peptides.build("?AAA?AA????A", aa)

mp = ProtoSyn.Forcefields.loadmap("amber03/aminoacids.yml");

# blist = ProtoSyn.genbonded(mol)

# top,graph,blist = ProtoSyn.Forcefields.assign(mol, ff, mp)#, blist)
top = ProtoSyn.Forcefields.assign(mol, ff, mp)#, blist)


# lr = LinkedResidue(source=aa["BKB"])
# push!(mol, lr)

# link = PS.bind(PS.Peptides.peptidebond, mol.residues[1], lr)
# link!==nothing && push!(mol, link)

# PS.update!(mol)




# dihedrals = PS.Peptides.finddihedrals(mol)
# PS.rotate!(state, mol, dihedrals[3], 0.5)


end
