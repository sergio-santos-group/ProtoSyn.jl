module Tst

include("../src/ProtoSyn.jl")

using .ProtoSyn
using .ProtoSyn.Builder

STEPWISE = false

if STEPWISE
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



    r13  = ProtoSyn.load("examples/r13.yml")
    r123 = ProtoSyn.load("examples/r123.yml")

    grammar = LGrammar{Char, String}()

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

    @pymol Builder.getvar(grammar, 'A')
    @pymol Builder.getvar(grammar, 'B')

    ProtoSyn.write(stdout, Builder.getvar(grammar, 'A'))
    ProtoSyn.write(stdout, Builder.getvar(grammar, 'B'))

    @show derivation = join(derive(grammar, "A", 4))

else

    using YAML
    yml = YAML.load(open("examples/data/protosugar.yml"))
    grammar = Builder.lgfactory(Float64, yml["protosugar"])

    @pymol Builder.getvar(grammar, "A")
    @pymol Builder.getvar(grammar, "B")

    @show derivation = derive(grammar, seq"A", 1)
    @show derivation = derive(grammar, derivation, 1)
    @show derivation = derive(grammar, derivation, 1)
    @show derivation = derive(grammar, derivation, 1)
    @show derivation = derive(grammar, derivation, 1)
    @show join(derivation)

end

# @pymol 

@pymol begin
    pose = build(grammar, derivation)
    sync!(pose)
end







end
