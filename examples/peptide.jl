module Tst

include("../src/ProtoSyn.jl")

using .ProtoSyn
using .ProtoSyn.Peptides
using .ProtoSyn.Builder

# grammar = Peptides.grammar()
grammar = Peptides.grammar(Float64)

# @pymol pose = build(grammar, "AAGASTASSE")
@pymol pose = build(grammar, seq"AAGASTASSE")

# println(pose.state)

# ProtoSyn.write(stdout, pose)

@pymol sync!(pose)
ProtoSyn.write(stdout, pose)

setss!(pose, SecondaryStructure[:linear])
@pymol sync!(pose)

# setss!(pose, SecondaryStructure[:helix])
# @pymol sync!(pose)



end