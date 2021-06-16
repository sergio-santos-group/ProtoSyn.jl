module Tst

include("../src/ProtoSyn.jl")

using YAML
using .ProtoSyn
using .ProtoSyn
using .ProtoSyn.Peptides
using .ProtoSyn.Calculators#.Forcefield
using .ProtoSyn.Drivers

grammar = Peptides.grammar(Float64)
pose = build(grammar, seq"AAAAAAAAAA")

setss!(pose, SecondaryStructure[:linear])
# @pymol sync!(pose)
sync!(pose)

ff = ForcefieldParameters("resources/Calculators/Forcefield/forcefield.yml")
resmaps = YAML.load_file("resources/Calculators/Forcefield/aminoacids.yml")

ftop = genff(pose.graph, ff, resmaps)
state = Calculators.State(pose.state)


const   FORCES = Val{true}
const NOFORCES = Val{false}


function myeval(s::Calculators.State, do_forces::Bool)
    if do_forces
        ProtoSyn.Calculators.eval!(s, ftop, FORCES)
    else
        ProtoSyn.Calculators.eval!(s, ftop, NOFORCES)
    end
    energy(state)    
end

sd = Drivers.SteepestDescent{Float64}(eval! = myeval)
sd.max_steps = 100_000

dstate = sd(state) do s,ds
    if ds.step%100 == 0
        println(ds.step,": ", Calculators.energy(s), " with stepsize=", ds.stepsize, " and maxf=", ds.max_force)
        copy!(pose.state, s)
        # @pymol pose
    end
end


end