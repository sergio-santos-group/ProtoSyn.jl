module Tst

include("../src/ProtoSyn.jl")

using .ProtoSyn
using .ProtoSyn.Peptides
using .ProtoSyn.Builder

# grammar = Peptides.grammar()
grammar = Peptides.grammar(Float64)

# @pymol pose = build(grammar, "AAGASTASSE")
pose = build(grammar, seq"AAQ")

# Try each one individually
selection = @resname "ALA"
# selection = @resname "ALA" | @resname "GLN"
# selection = @resname "ALA" & @resname "GLN"
# selection = (@resname "ALA" | @resname "GLN") & @resname "ALA"

println(selection(pose.graph))
exit(1)

sync!(pose)

selection = @resname "ALA"
# selection = @resname "ALA" | @resname "GLN"
# selection = @resname "ALA" & @resname "GLN"
# selection = (@resname "ALA" | @resname "GLN") & @resname "ALA"

println(selection(pose.graph))

exit(1)
# println(pose.state)

# ProtoSyn.write(stdout, pose)

ProtoSyn.write(stdout, pose)

setss!(pose, SecondaryStructure[:linear])
@pymol sync!(pose)

# setss!(pose, SecondaryStructure[:helix])
# @pymol sync!(pose)



end