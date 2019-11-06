at1 = AtomMetadata(index=1, name="atom1", symbol="W")
at2 = AtomMetadata(index=2, name="atom2", symbol="X")
at3 = AtomMetadata(index=3, name="atom3", symbol="Y")
at4 = AtomMetadata(index=4, name="atom4", symbol="Z")
at5 = AtomMetadata(index=5, name="atom5", symbol="A")
at6 = AtomMetadata(index=6, name="atom6", symbol="A")
at7 = AtomMetadata(index=7, name="atom7", symbol="A")
at8 = AtomMetadata(index=8, name="atom8", symbol="A")

r1 = ResidueMetadata(index=1, name="res1")
r2 = ResidueMetadata(index=2, name="res2")
r3 = ResidueMetadata(index=3, name="res3")

m1 = MoleculeMetadata(index=1, name="mol1")

push!(r1, at1, at2, at3)
push!(r2, at4, at5)
push!(r3, at6, at7, at8)

m1.children = [r1,r2,r3]

println(m1)
println(r1)
println(r2)
println(r3)




println(findfirst(r->r===r2, m1.children))

println(delete!(m1, r1))
# println(insert!(m1, 3, r1).children)

# iterator = IndicesIterator([at1,at2,at3], [r1, r2])
# for i in iterator
#     println(">> ", i.index)
# end

#x,m = readpdb("ProtoSyn.jl/private/2_peptides.pdb", true)
#x,m = readpdb("../../private/2_peptides.pdb", true)
#println(size(x))
# println(m[:atoms])
# println(m[:residues])
#println(m[:molecules])
# println(isa(m, Nothing))

