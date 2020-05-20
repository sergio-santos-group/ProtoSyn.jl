module LSystem

using ..ProtoSyn

export StochasticRule, LGrammar
export derive, build

struct StochasticRule{K,V}
    p::Float64
    # rule::T
    source::K
    production::V
    StochasticRule(p::Float64, rule::Pair{K,V}) where {K,V} = begin
        (p > 1 || p < 0) && error("Invalid probability")
        new{K,V}(p, rule.first, rule.second)
    end
end


struct LGrammar{K,V}
    rules::Dict{K, Vector{StochasticRule{K,V}}}
    variables::Dict{K, Fragment}
    operators::Dict{K, Function}
end


LGrammar{K,V}() where {K,V} = LGrammar(
    Dict{K, Vector{StochasticRule{K,V}}}(),
    Dict{K, Fragment}(),
    Dict{K, Function}()
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
    ruleset[i].production
end


getvar(g::LGrammar, k) = begin
    !haskey(g.variables, k) && throw(KeyError(k))
    g.variables[k]
end


getop(g::LGrammar, k) = begin
    !haskey(g.operators, k) && throw(KeyError(k))
    g.operators[k]
end


Base.push!(g::LGrammar{K,V}, r::StochasticRule{K,V}) where {K,V} = begin
    push!(
        get!(g.rules, r.source, []),
        r
    )
    g
end


Base.setindex!(g::LGrammar{K,V}, x::Fragment, key::K) where {K,V} = begin
    g.variables[key] = x
    g
end


Base.setindex!(g::LGrammar{K,V}, x::Function, key::K) where {K,V} = begin
    g.operators[key] = x
    g
end


derive(lg::LGrammar, axiom) = [
    x
    for r in axiom
        for x in getrule(lg,r) ]

derive(lg::LGrammar, axiom, niter::Int) = begin
   niter > 0 ? derive(lg, derive(lg, axiom), niter-1) : axiom
end


isop(lg::LGrammar, k) = haskey(lg.operators, k)
isvar(lg::LGrammar, k) = haskey(lg.variables, k)

function fragment(::Type{T}, grammar::LGrammar, derivation) where {T<:AbstractFloat}

    state = State{T}()
    seg = Segment("derivation", -1)

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
            
            frag2 = copy(frag)

            push!(seg, frag2.graph.items...)
            append!(state, frag2.state)
            if parent !== nothing
                join = pop!(opstack)
                #println(join)
                join(parent, frag2)
            end
            parent = frag2
        end
    end
    reindex(seg)
    seg.id = state.id = ProtoSyn.genid()

    return Pose(seg, state)
end

build(grammar::LGrammar, derivation) = build(Float64, grammar, derivation)

build(::Type{T}, grammar::LGrammar, derivation) where {T<:AbstractFloat} = begin
    top = Topology("unknown", 1)
    state = State{T}()
    state.id = top.id
    pose = Pose(top, state)

    if !isempty(derivation)
        frag = fragment(T, grammar, derivation)
        append(pose, frag)
        # reindex(top)
        ProtoSyn.request_i2c(state; all=true)
    end
    pose
end

end