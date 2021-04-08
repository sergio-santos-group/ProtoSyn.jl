module Tst

include("../src/ProtoSyn.jl")

using YAML
using .ProtoSyn
using .ProtoSyn.Builder
using .ProtoSyn.Peptides
using .ProtoSyn.Drivers
using .ProtoSyn.Calculators

grammar = Peptides.grammar(Float64)
pose = build(grammar, seq"AAAAAAAAAA")

setss!(pose, SecondaryStructure[:linear])
@pymol sync!(pose)

ff = ForcefieldParameters("resources/Calculators/Forcefield/forcefield.yml")
resmaps = YAML.load_file("resources/Calculators/Forcefield/aminoacids.yml")

ftop = genff(pose.graph, ff, resmaps)
state = pose.state


const   FORCES = Val{true}
const NOFORCES = Val{false}


function myeval(s::InternalRepresentation, ::Bool)
    ProtoSyn.request_i2c!(s.state)
    ProtoSyn.sync!(s.state, pose.graph)

    st = s.state
    x = st.x
    for i=1:st.size
        t = st[i].t
        x[i,1] = t[1]
        x[i,2] = t[2]
        x[i,3] = t[3]
    end
    ProtoSyn.Calculators.eval!(s.state, ftop, NOFORCES)
    energy(state)    
end


dihedrals = [
    r[dihd]
    for r in eachresidue(pose.graph)
    for dihd in (
        ProtoSyn.Peptides.DihedralTypes.phi,
        ProtoSyn.Peptides.DihedralTypes.psi
    )
]

function mysampler(s::InternalRepresentation)
    for dihedral in dihedrals
        if rand() < 0.1
            setdihedral!(s.state, dihedral, randn())
        end
    end
    #ProtoSyn.request_i2c!(pose.state)
    #sync!(pose)
end




sd = Drivers.MonteCarlo{Float64}(
    eval! = myeval,
    sample! = mysampler,
    temperature = 100.0
)
sd.max_steps = 100
# sd.max_steps = 10

dstate = sd(InternalRepresentation(state)) do s,ds
    if ds.step%100 == 0
        println(ds.step,": ", energy(s))
        #copy!(pose.state, s)
        @pymol pose
    end
end


end