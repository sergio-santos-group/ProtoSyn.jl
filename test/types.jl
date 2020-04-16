
lib = ProtoSyn.Peptides.load(ProtoSyn.Peptides.resource_dir)


# mol, state = ProtoSyn.Peptides.build("AAE", lib)

# foreach(r->println("$(r.name):$(r.offset)-$(length(r))"), mol.residues)

# open("test.pdb", "w") do io
#     write(io, mol, state, ProtoSyn.PDB)

#     replace!(mol, state, mol.residues[2]=>lib["SER"])

#     write(io, mol, state, ProtoSyn.PDB)
#     foreach(r->println("$(r.name):$(r.offset)-$(length(r))"), mol.residues)
#     println("state size=$(state.size), capacity=$(state.capacity)")
# end


#mol, state = ProtoSyn.Peptides.build("EFDVILKAAGANKVAVIKAV", lib)
#mol, state = ProtoSyn.Peptides.build("G"^20, lib)
#println(join(keys(ProtoSyn.Peptides.one_2_three), ""))
mol, state = ProtoSyn.Peptides.build("MKPQIHEWSTC?DALYVRGFN", lib)


#exit(0)
# r = ProtoSyn.rotmat(rand(3), rand())
# ProtoSyn.rotate!(state.coords, r)
# state.coords .+= rand()*5


pymol = ProtoSyn.XMLRPC.ServerProxy("http://localhost", 9123)
pymol.delete("all")
io = IOBuffer()

write(io, mol, state, ProtoSyn.PDB)
pymol.read_pdbstr(String(take!(io)), "helix_ref")

ϕ, ψ = deg2rad(-60.0), deg2rad(-45.0)
ProtoSyn.set!(state, mol, 1:20, ("-C","N","CA","C"), ϕ)
ProtoSyn.set!(state, mol, 1:20, ("N","CA","C","+N"), ψ)

write(io, mol, state, ProtoSyn.PDB)
pymol.read_pdbstr(String(take!(io)), "helix_ref")

println("BEFORE EVERYTHING ", state.size)
foreach(println, eachcol(state.coords))

ProtoSyn.Peptides.replace!(mol, state, mol.residues[10]=>lib["SER"])
ProtoSyn.Peptides.replace!(mol, state, mol.residues[15]=>lib["PRO"])
write(io, mol, state, ProtoSyn.PDB)
pymol.read_pdbstr(String(take!(io)), "helix")

ProtoSyn.Peptides.replace!(mol, state, mol.residues[15]=>lib["GLY"])
write(io, mol, state, ProtoSyn.PDB)
pymol.read_pdbstr(String(take!(io)), "helix2")


exit(1)

# open("helix_ref.pdb", "w") do io
#     # helix
#     ϕ, ψ = deg2rad(-60.0), deg2rad(-45.0)
#     ProtoSyn.set!(state, mol, 1:20, ("-C","N","CA","C"), ϕ)
#     ProtoSyn.set!(state, mol, 1:20, ("N","CA","C","+N"), ψ)
#     write(io, mol, state, ProtoSyn.PDB)
# end


# open("helix.pdb", "w") do io
#     replace!(mol, state, mol.residues[10]=>lib["SER"])
#     replace!(mol, state, mol.residues[15]=>lib["PRO"])
#     write(io, mol, state, ProtoSyn.PDB)
# end
# open("
# .pdb", "w") do io
#     replace!(mol, state, mol.residues[15]=>lib["GLY"])
#     write(io, mol, state, ProtoSyn.PDB)
# end

# write(io, mol, state, ProtoSyn.PDB)
# s = String(take!(io))
# pymol.read_pdbstr(s, "test")

# mol, state = ProtoSyn.Peptides.build("G"^5, lib)

dihedrals = filter(
    d->d.type != ProtoSyn.Peptides.DihedralTypes.omega,
    ProtoSyn.Peptides.finddihedrals(mol; sidechain=false)
)

foreach(println, dihedrals)

foreach(println, mol.links)
println(mol.bonds)

write(io, mol, state, ProtoSyn.PDB)
println(String(take!(io)))
# println(state.coords')
foreach(println, eachcol(state.coords))

foreach(println, mol.residues)
foreach(println, ProtoSyn.iterbyatom(mol.))

# for i in eachcol(state.coords)
#     println(i)
#     # println(state.coords[:,i])
# end
#for atom in ProtoSyn.iterbyatom(mol)
exit(1)


function sampler(s::ProtoSyn.State)
    for dihedral in dihedrals
        (rand() < 0.1) && ProtoSyn.rotate!(s,mol,dihedral, 0.1*rand())
        #ts = mol.traverse_state
        #println(ts.size, dihedral)
        #println(ts.indices[1:ts.size])
        #ProtoSyn.update!(mol.traverse_state, dihedral, mol.bonds)
    end
end

driver = ProtoSyn.Drivers.MonteCarlo(
    (s::ProtoSyn.State, ::Bool) -> 0.0,
    sampler,
    1000,
    1.0
)
write(io, mol, state, ProtoSyn.PDB)
pymol.read_pdbstr(String(take!(io)), "mcarlo")
# for i = 1:10
#     sampler(state)
#     write(io, mol, state, ProtoSyn.PDB)
#     pymol.read_pdbstr(String(take!(io)), "mcarlo")
# end
pymol.do("axes on")
driver(state) do st,ds
    if ds.step%50 == 0
        write(io, mol, state, ProtoSyn.PDB)
        #s = String(take!(io))
        pymol.read_pdbstr(String(take!(io)), "mcarlo")
    end
end