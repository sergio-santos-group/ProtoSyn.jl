module Tst

include("../src/ProtoSyn.jl")

using .ProtoSyn
using .ProtoSyn.Sugars
using .ProtoSyn.Builder

grammar = Sugars.grammar(Float64, "amylose")

@pymol pose = build(grammar, seq"AAAAAAAAAA")
@pymol sync!(pose)

ProtoSyn.write(stdout, pose)

end