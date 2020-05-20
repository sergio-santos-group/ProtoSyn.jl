module Tst

include("../src/ProtoSyn.jl")

using .ProtoSyn
using .ProtoSyn.LSystem

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



r13  = read(Float64, "examples/r13.yml",  ProtoSyn.YML)
r123 = read(Float64, "examples/r123.yml", ProtoSyn.YML)

grammar = LGrammar{Char,String}()

# add variables
grammar['A'] = ProtoSyn.fragment(r13)
grammar['B'] = ProtoSyn.fragment(r123)

# add operators
grammar['α'] = α31
grammar['β'] = α21

# add rules
# push!(grammar, StochasticRule(1.00, 'B' => "AαA"))
push!(grammar, StochasticRule(0.75, 'A' => "AαA"))
push!(grammar, StochasticRule(0.25, 'A' => "B[αA]βA"))



mytoolbelt = ReactionToolbelt(
    ()->(),
    ()->(),
    (s)->s[1,"C1"],
)

@pymol LSystem.getvar(grammar, 'A')
@pymol LSystem.getvar(grammar, 'B')

write(stdout, LSystem.getvar(grammar, 'A'))
write(stdout, LSystem.getvar(grammar, 'B'))

@show derivation = derive(grammar, "A", 3)

@pymol begin
    # @show derivation = "B[αAαAαAαA]βB[αA]βAαAαA"
    # top = Topology("UNK", 1)
    # state = State{Float64}()
    # state.id = top.id
    # pose = Pose(top, state)

    # frag = build(grammar, derivation)
    # frag.graph.name = "A"
    # frag.graph.id = 1
    # append(pose, frag, mytoolbelt)
    # reindex(top)
    # ProtoSyn.request_i2c(pose.state; all=true)
    pose = build(grammar, derivation)
    sync!(pose)
end

end
