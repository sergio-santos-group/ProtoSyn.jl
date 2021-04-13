module Tst

include("../src/ProtoSyn.jl")

using YAML
using .ProtoSyn
using .ProtoSyn
using .ProtoSyn.Peptides
using .ProtoSyn.Calculators#.Forcefield

grammar = Peptides.grammar(Float64)
pose = build(grammar, seq"AAAAAAAAAA")

setss!(pose, SecondaryStructure[:linear])
@pymol sync!(pose)

ff = ForcefieldParameters("resources/Calculators/Forcefield/forcefield.yml")
resmaps = YAML.load_file("resources/Calculators/Forcefield/aminoacids.yml")
# blist = genbonded(pose.graph,3)


ftop = genff(pose.graph, ff, resmaps)
state = Calculators.State(pose.state)

const   FORCES = Val{true}
const NOFORCES = Val{false}

@noinline function f1(n::Int, dof)
    for i=1:n
        # ProtoSyn.Calculators.eval!(state, ftop, dof)
        eval!(state, ftop, dof)
    end
    println(energy(state))
end

@time f1(    1, NOFORCES)
@time f1(    1,   FORCES)
@time f1(1000, NOFORCES)
@time f1(1000,   FORCES)

function sd(n)
    λ = 0.00005
    eOld = Inf64

    for i=1:n
        ProtoSyn.Calculators.eval!(state, ftop, FORCES)
        eNew = energy(state)
        copy!(pose.state, state)
        
        if i%1==0
            println("Step ", i, " E = ", eNew, " with λ=", λ)
            for (k,v) in state.e
                # println("   ", name(k),"{", k.parameters[2],"} = ", v)
                println("   ", k, " = ", v)
            end
            # @pymol pose
        end
        @. state.x += λ*state.f
        if eNew < eOld
            λ*=1.05
        else
            λ*=0.5
        end
        eOld = eNew

    end
end
sd(10)

end