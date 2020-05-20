module LSystem

using ..ProtoSyn

# include("graph.jl")
# include("../ProtoSyn.jl")
# using ..ProtoSyn

mutable struct LNode <: ProtoSyn.AbstractDigraph
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
    variables::Dict{Char,Fragment}
    operators::Dict{Char, Function}
end

LGrammar() = LGrammar(
    Dict{Char, StochasticRule}(),
    Dict{Char, Fragment}(),
    Dict{Char, Function}()
)


getrule(g::LGrammar, r) = begin
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
    ruleset[i].rule
end

getvar(g::LGrammar, k) = begin
    !haskey(g.variables, k) && throw(KeyError(k))
    g.variables[k]
end

getop(g::LGrammar, k) = begin
    !haskey(g.operators, k) && throw(KeyError(k))
    g.operators[k]
end

addrule!(g::LGrammar, r::Pair{Char, StochasticRule}) = begin
    rules = get!(g.rules, r.first, ())
    g.rules[r.first] = (rules..., r.second)
    g
end

derive(lg::LGrammar, axiom::String) = join(getrule(lg,c) for c in axiom)
derive(lg::LGrammar, axiom::String, niter::Int) = begin
   niter > 0 ? derive(lg, derive(lg, axiom), niter-1) : axiom
end

function build(derivation::String)
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

Base.getindex(lg::LGrammar, i) = begin
    haskey(lg.operators, i) ? getop(lg, i) : getvar(lg, i)
end
isop(lg::LGrammar, k) = haskey(lg.operators, k)
isvar(lg::LGrammar, k) = haskey(lg.variables, k)

function build(grammar::LGrammar, derivation::String)

    seg = Segment(derivation, -1)
    state = State{Float64}()

    
    stack = Fragment[]
    opstack = Function[]
    parent::Opt{Fragment} = nothing

    for letter in derivation
        if isop(grammar, letter)
            op = getop(grammar, letter)
            push!(opstack, op)
        elseif letter == '['
            push!(stack, parent)
        elseif letter == ']'
            parent = pop!(stack)
        elseif isvar(grammar, letter)
            frag = getvar(grammar, letter)
            # join = isempty(opstack) ? defop : pop!(opstack)
            

            #r = copy(frag.graph)
            #s = copy(frag.state)
            frag2 = copy(frag)

            push!(seg, frag2.graph.items...)
            append!(state, frag2.state)
            if parent !== nothing
                join = pop!(opstack)
                println(join)
                join(parent, frag2)
            end
            parent = frag2
        end
    end
    reindex(seg)
    seg.id = state.id = ProtoSyn.genid()

    return Pose(seg, state)
end


# grammar = LGrammar()
# addrule!(grammar, '1' => StochasticRule(1.0, "11"))
# addrule!(grammar, '0' => StochasticRule(1.5, "1[0]0"))
# addrule!(grammar, '0' => StochasticRule(0.5, "00"))

# for i=1:4
#     @show derive(grammar, "0", i)
# end

# sugar = LGrammar()
# addrule!(sugar, 'B' => StochasticRule(1.00, "AαA"))
# addrule!(sugar, 'A' => StochasticRule(0.75, "AαA"))
# addrule!(sugar, 'A' => StochasticRule(0.25, "B[βA]αA"))

# for i=1:4
#     @show derive(sugar, "B", i)
# end

#--------------------

# struct LLGrammar{T}
#     rules::Dict{T,Tuple{Vararg{StochasticRule}}}
#     variables::Dict{T,Fragment}
#     operators::Dict{T, Function}
# end

# LLGrammar{T}() where T = LLGrammar{T}(
#     Dict{T,StochasticRule}(),
#     Dict{T,Fragment}(),
#     Dict{T, Function}()
# )
# LLGrammar() = LLGrammar{Char}()

# Base.push!(g::LLGrammar, r::Pair{T, StochasticRule}) where T = begin
#     rules = get!(g.rules, r.first, ())
#     g.rules[r.first] = (rules..., r.second)
#     g
# end

# Base.get(g::LLGrammar, r) = begin
#     if !haskey(g.rules, r)
#         return r
#     end
    
#     i = 1
#     p = rand()
#     ruleset = g.rules[r]
#     x = ruleset[1].p
#     while x < p
#         x += ruleset[(i+=1)].p
#     end
#     ruleset[i].rule
# end

# derive(lg::LLGrammar{Char}, axiom::String) = join(get(lg,c) for c in axiom)
# derive(lg::LLGrammar{Char}, axiom::String, niter::Int) = begin
#    niter > 0 ? derive(lg, derive(lg, axiom), niter-1) : axiom
# end


end