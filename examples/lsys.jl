module Tst

include("../src/ProtoSyn.jl")

using .ProtoSyn


function α31(f1::Fragment, f2::Fragment)
    r1 = f1.graph[end]
    r2 = f2.graph[1]

    if hasparent(r2)
        error("r2 is already connected")
    else
        c3 = r1["C3"]
        c1 = r2["C1"]
        bond(c3,c1)
        setparent!(c1, c3)  # c3->c1
        setparent!(r2, r1)  # r1->r2
    end
end
# α31(s1::Segment, s2::Segment) = α31(s1[end], s2[1])
# α31(r1::Residue, s2::Segment) = α31(r1, s2[1])
# α31(s1::Segment, r2::Residue) = α31(s1[end], r2)


function α21(f1::Fragment, f2::Fragment)
    r1 = f1.graph[end]
    r2 = f2.graph[1]

    if hasparent(r2)
        error("r2 is already connected")
    else
        c2 = r1["C2"]
        c1 = r2["C1"]
        bond(c2,c1)
        setparent!(c1, c2)  # c2->c1
        setparent!(r2, r1)  # r1->r2

        state = f2.state
        st = state[r1["C1"]]
        st.b = 2.5
        st.θ = deg2rad(135)
        st.ϕ = deg2rad(0)
    end
end
# α21(s1::Segment, s2::Segment) = α21(s1[end], s2[1])
# α21(r1::Residue, s2::Segment) = α21(r1, s2[1])
# α21(s1::Segment, r2::Residue) = α21(s1[end], r2)



r13  = read(Float64, "examples/r13.yml", ProtoSyn.YML)
r123 = read(Float64, "examples/r123.yml", ProtoSyn.YML)

grammar = ProtoSyn.LSystem.LGrammar()

# add variables
grammar.variables['A'] = ProtoSyn.fragment(r13)
grammar.variables['B'] = ProtoSyn.fragment(r123)


# add rules
# ProtoSyn.LSystem.addrule!(grammar, 'B' => ProtoSyn.LSystem.StochasticRule(1.00, "AαA"))
ProtoSyn.LSystem.addrule!(grammar, 'A' => ProtoSyn.LSystem.StochasticRule(0.75, "AαA"))
ProtoSyn.LSystem.addrule!(grammar, 'A' => ProtoSyn.LSystem.StochasticRule(0.25, "B[αA]βA"))


# add ...
grammar.operators['α'] = α31
grammar.operators['β'] = α21


mytoolbelt = ReactionToolbelt(
    ()->(),
    ()->(),
    (s)->s[1,"C1"],
)

@pymol grammar.variables['A']
@pymol grammar.variables['B']

write(stdout, grammar.variables['A'])
write(stdout, grammar.variables['B'])

@show derivation = ProtoSyn.LSystem.derive(grammar, "A", 3)
@pymol begin
    @show derivation = "B[αAαAαAαA]βB[αA]βAαAαA"
    top = Topology(derivation, 1)
    state = State{Float64}()
    state.id = top.id
    pose = Pose(top, state)

    frag = ProtoSyn.LSystem.build(grammar, derivation)
    frag.graph.name = "A"
    frag.graph.id = 1
    append(pose, frag, mytoolbelt)
    reindex(top)
    ProtoSyn.request_i2c(pose.state; all=true)
    sync!(pose)
end

end
