# module Tst

include("../src/ProtoSyn.jl")

using ProtoSyn
using ProtoSyn.Peptides
using ProtoSyn.Builder

# grammar = Peptides.grammar()
grammar = Peptides.grammar(Float64)
println("Loaded all grammar.")

# @pymol pose = build(grammar, "AAGASTASSE")
pose = build(grammar, seq"AAQG")
println("Done building.")

setss!(pose, SecondaryStructure[:linear])
sync!(pose)

io = open("./teste.pdb", "w")
print_selection(io, pose, (@resname r"GL." & @atomsymb "C")(pose.graph))
close(io)
println("Done printing selection to file.")

# println(selection(pose.graph))
# exit(1)

# sync!(pose)

# # println(pose.state)

# # ProtoSyn.write(stdout, pose)

# ProtoSyn.write(stdout, pose)

# setss!(pose, SecondaryStructure[:linear])
# @pymol sync!(pose)

# setss!(pose, SecondaryStructure[:helix])
# @pymol sync!(pose)



# end