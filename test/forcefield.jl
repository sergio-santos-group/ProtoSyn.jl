module Tst

include("../src/Calculators/Forcefield/Forcefield.jl")
using .Forcefield

ff = loadparm("resources/Calculators/Forcefield/protoff.yml")

end


# module Tst

# #include("../src/ProtoSyn.jl")

# # using .ProtoSyn
# # using .ProtoSyn.Sugars
# # using .ProtoSyn.Builder

# # grammar = Sugars.grammar(Float64, "amylose")

# # # @pymol pose = build(grammar, seq"AAAAAAAAAA")
# # # ProtoSyn.write(stdout, pose)

# #using .ProtoSyn.Calculators.Forcefield

# # pose = build(grammar, seq"AA")
# # @pymol sync!(pose)

# # blist = genbonded(pose.graph[1], 3)
# include("../src/Calculators/Forcefield/Forcefield.jl")
# using .Forcefield
# ff = loadparm("resources/Calculators/Forcefield/protoff.yml")

# end