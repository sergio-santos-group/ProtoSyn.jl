module LSystem

include("graph.jl")

mutable struct LNode <: AbstractDigraph
    parent::Union{LNode,Nothing}
    children::Vector{LNode}
    visited::Bool
    name::Char
    index::Int
    LNode(p::Union{LNode,Nothing}, name::Char, index::Int) = begin
        node = new(nothing, LNode[], false, name, index)
        if p !== nothing
            node.parent = p
            push!(p.children, node)
        end
        node
    end
end

struct StochasticRule
    p::Float64
    rule::String
end


struct LGrammar
    rules::Dict{Char,Tuple{Vararg{StochasticRule}}}
    # LGrammar(rules::Pair{String,StochasticRule}...)
end
LGrammar() = LGrammar(Dict{Char,StochasticRule}())

Base.get(g::LGrammar, r) = begin
    if !haskey(g.rules, r)
        return r
    end
    
    i = 1
    p = rand()
    ruleset = g.rules[r]
    x = ruleset[1].p
    while x < p
        x += ruleset[(i+=1)].p
    end
    g.rules[r][i].rule
end

Base.push!(g::LGrammar, r::Pair{Char, StochasticRule}) = begin
    rules = get!(g.rules, r.first, ())
    g.rules[r.first] = (rules..., r.second)
    g
end

(lg::LGrammar)(axiom::String, niter::Int) = begin
    derivations = String[axiom]
    for i=1:niter
        axiom = derive(lg,axiom)
        push!(derivations, axiom)
    end
    derivations
end

derive(lg::LGrammar, axiom::String) = join([get(lg,c) for c in axiom])

function embed(derivation::String)
    graph = LNode[]
    stack = LNode[]
    node = nothing
    for letter in derivation
        if letter == '['
            push!(stack, node)
        elseif letter == ']'
            node = pop!(stack)
        else
            node = LNode(node, letter, length(graph)+1)
            push!(graph, node)
        end
    end
    graph
end



grammar = LGrammar()
push!(grammar, '1' => StochasticRule(1.0, "11"))
push!(grammar, '0' => StochasticRule(0.5, "1[0]0"))
push!(grammar, '0' => StochasticRule(0.5, "00"))

# grammar = LGrammar(
#     Dict(
#         '1' => (1.0, "11"),
#         '0' => ((0.5, "1[0]0"), (0.5, "0"))
#     )
# )

@show grammar("0", 4)

end