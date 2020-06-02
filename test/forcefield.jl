# module Tst

# include("../src/Calculators/Forcefield/Forcefield.jl")
# using .Forcefield

# ff = loadparm("resources/Calculators/Forcefield/protoff.yml")

# end


module Tst

include("../src/ProtoSyn.jl")

using .ProtoSyn
using .ProtoSyn.Builder
using .ProtoSyn.Peptides
using .ProtoSyn.Calculators.Forcefield

grammar = Peptides.grammar(Float64)
pose = build(grammar, seq"AA")

# setss!(pose, SecondaryStructure[:linear])
# @pymol sync!(pose)

ff = ForcefieldParameters("resources/Calculators/Forcefield/protoff.yml")
blist = genbonded(pose.graph,3)

using YAML
resmaps = YAML.load_file("resources/Calculators/Forcefield/aminoacids.yml")
# atypes = assigntypes(pose.graph[1], resmaps)

ftop = genff(pose.graph[1], ff, resmaps)
state = Forcefield.State(pose.state)

const   FORCES = Val{true}
const NOFORCES = Val{false}

# ProtoSyn.Calculators.eval!(
#     state,
#     ftop.components[Forcefield.HarmonicPotential{Int64,2,Float64}],
#     FORCES
#     )
# ProtoSyn.Calculators.eval!(
#     state,
#     ftop.components[Forcefield.HarmonicPotential{Int64,2,Float64}],
#     NOFORCES
#     )

@noinline function f1(n::Int, dof)
    # e = 0.0
    # for i=1:n
    #     e += i
    # end
    bonds = ftop.components[Forcefield.HarmonicPotential{Int64,2,Float64}]
    #println(typeof(bonds))
    for i=1:n
        #e += 
        ProtoSyn.Calculators.eval!(state, bonds, dof)
    end
    
end

# f1(10,   FORCES)
# f1(10, NOFORCES)
# f1(100000,   FORCES)
# f1(100000, NOFORCES)
# f(n,f) = @time f1(n,f)

@time f1(    1,   FORCES)
@time f1(    1, NOFORCES)
@time f1(100000,   FORCES)
@time f1(100000, NOFORCES)

end