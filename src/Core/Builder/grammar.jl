export StochasticRule, LGrammar
export derive, build, lgfactory, fromfile
export @seq_str

struct StochasticRule{K, V}
    p::Float64
    source::K
    production::V
    StochasticRule(p::Float64, rule::Pair{K,V}) where {K,V} = begin
        (p > 1 || p < 0) && error("Invalid probability")
        new{K,V}(p, rule.first, rule.second)
    end
end
Base.show(io::IO, sr::StochasticRule) = begin
    println(io, sr.source, "(p=", sr.p, ") -> ", sr.production)
end

mutable struct LGrammar{T <: AbstractFloat, K, V}
    rules::Dict{K, Vector{StochasticRule{K,V}}}
    variables::Dict{K, Fragment}
    operators::Dict{K, Function}
    defop::Opt{Function}
end


LGrammar{T, K, V}() where {T <: AbstractFloat, K, V} = LGrammar{T, K, V}(
    Dict{K, Vector{StochasticRule{K,V}}}(),
    Dict{K, Fragment}(),
    Dict{K, Function}(),
    nothing
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

hasrule(g::LGrammar, r) = haskey(g.rules, r)

Base.push!(g::LGrammar{T, K, V}, r::StochasticRule{K,V}) where {T <: AbstractFloat, K, V} = begin
    push!(
        get!(g.rules, r.source, []),
        r
    )
    g
end


Base.setindex!(g::LGrammar{T, K, V}, x::Fragment, key::K) where {T <: AbstractFloat, K, V} = begin
    g.variables[key] = x
    g
end


Base.setindex!(g::LGrammar{T, K, V}, x::Function, key::K) where {T <: AbstractFloat, K, V} = begin
    g.operators[key] = x
    g
end


derive(lg::LGrammar, axiom) = begin
    derivation = []
    for x in axiom
        if hasrule(lg, x)
            append!(derivation, getrule(lg, x))
        else
            push!(derivation, x)
        end
    end
    derivation
end
# [
#     x
#     for r in axiom
#         for x in getrule(lg,r) ]

derive(lg::LGrammar, axiom, niter::Int) = begin
   niter > 0 ? derive(lg, derive(lg, axiom), niter-1) : axiom
end


isop(lg::LGrammar, k) = haskey(lg.operators, k)
isvar(lg::LGrammar, k) = haskey(lg.variables, k)


#region fragment ----------------------------------------------------------------

export isfragment

isfragment(p::Pose) = !(hascontainer(p.graph) || isempty(p.graph))


"""
    fragment(grammar::LGrammar{T, K, V}, derivation) where {T <: AbstractFloat, K, V}
    
Create and return a new fragment (`Pose` instance with just a single `Segment`)
using the given `derivation` sequence on the provided `grammar` instructions.
The main purpose of fragments is to be temporary carriers of information, such
as during the building process of a new peptide from a sequence. Therefore,
these structures often don't have any real meaning and are, as such, deprived of
a root/origin for the graph. Actual structures should instead be of the slightly
more complete type `Pose`.

!!! note
    A fragment does not contain a `Topology` instance.

# Examples
```jldoctest
julia> frag = fragment(res_lib, seq"AAA")
```
"""
function fragment(grammar::LGrammar{T, K, V}, derivation) where {T <: AbstractFloat, K, V}

    state = State(T)
    seg = Segment("UNK", 1)
    seg.code = 'A'

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

            push!(seg, frag2.graph.items...) # Appending the residues to the segment
            append!(state, frag2.state)      # Merging the 2 states
            if parent !== nothing
                join = isempty(opstack) ? grammar.defop : pop!(opstack)
                join(parent.graph[end], frag2) # Adding ascendents and bonds correctly
            end
            parent = frag2
        end
    end
    reindex(seg)
    seg.id = state.id = ProtoSyn.genid()

    return Pose(seg, state)
end

# fragment(grammar::LGrammar, derivation) = fragment(Units.defaultFloat, grammar, derivation)

#endregion fragment


"""
    opfactory(args) -> Function

Creates the operation function (as closures) given the input arguments `args`
(normally read from a grammar file). 

Note: The operation functions are responsible for setting the internal
coordinates of residues in the system when connecting, building and manipulating
poses.
"""
function opfactory(args)
    # Note: residue_index is the index on the fragment/pose ('f2') used to
    # actually connect to 'r1'
    return function(r1::Residue, f2::Union{Fragment, Pose}; residue_index = 1)
        # residue_index is of f2
        if typeof(f2) == Fragment
            r2 = f2.graph[residue_index]
        else
            r2 = f2.graph[1][residue_index]
        end
        ProtoSyn.join(r1, args["residue1"], r2, args["residue2"]) # Connects specifically C to N (in case of proteins) -> Adds bonds and sets parents
        state = f2.state
        
        if haskey(args, "presets")
            for (atname, presets) in args["presets"]
                !(atname in r2) && continue # Proline has no H
                atomstate = state[r2[atname]]
                for (key, value) in presets
                    setproperty!(atomstate, Symbol(key), value)
                end
            end
        end
        
        if haskey(args, "offsets")
            for (atname, offset) in args["offsets"]
                setoffset!(state, r2[atname], offset)
            end
        end
    end
end



"""
    lgfactory([T=Float64,] template::Dict) -> LGrammar{T, String, Vector{String}}

Create an `LGrammar` instance from the contencts of a `template` Dict (normally
read from a grammar file).

# Example
```yml
amylose:
  rules:
    A:
      - {p: 0.75, production: [A,α,A]}
      - {p: 0.25, production: [B,"[",α,A,"]",β,A]}
  variables:
    A: resources/Sugars/GLC14.pdb
  operators:
    α14:
      residue1: C1
      residue2: O4
      presets:
        O4: {θ: 127, ϕ: -90.712, b: 1.43}
        C4: {θ: 100, ϕ: -103.475}
      offsets:
        H4: 0
  defop: α14
```
"""
function lgfactory(::Type{T}, template::Dict) where T
    grammar = LGrammar{T, String, Vector{String}}()

    vars = template["variables"]
    for (key, name) in vars
        filename = joinpath(ProtoSyn.resource_dir, name)
        @info "Loading variable '$key' from $filename"
        pose = ProtoSyn.load(T, filename)
        grammar[key] = fragment(pose)
    end

    ops = template["operators"]
    for (opname, opargs) in ops
        if haskey(opargs, "presets")
            for presets in values(opargs["presets"])
                for (k, v) in presets
                    presets[k] = ProtoSyn.Units.tonumber(T, v)
                end
            end
        end

        if haskey(opargs, "offsets")
            offsets = opargs["offsets"]
            for (k,v) in offsets
                offsets[k] = ProtoSyn.Units.tonumber(T, v)
            end
        end

        @info "Loading operator $opname"
        grammar[opname] = opfactory(opargs)
    end
    grammar.defop = getop(grammar, template["defop"])

    if haskey(template, "rules")
        for (key,rules) in template["rules"]
            @info "Loading productions for rule $key"
            for rule in rules
                sr = StochasticRule(rule["p"], key => rule["production"])
                push!(grammar, sr)
                @info "  $sr"
            end
        end
    end
    grammar
end

lgfactory(template::Dict) = lgfactory(Float64, template)


"""
    fromfile(::Type{T}, filename::AbstractString, key::String) -> LGrammar{String,Vector{String}}

Create an `LGrammar` instance from the contencts of a grammar file (in YML
format).

# Examples
```jldoctest
julia> lgrammar = fromfile(Float64, filename, "peptide")

julia> lgrammar = fromfile(filename, "peptide")
```
"""
function fromfile(::Type{T}, filename::AbstractString, key::String) where T
    open(filename) do io
        @info "loading grammar from file $filename"
        yml = YAML.load(io)
        lgfactory(T, yml[key])
    end
end

fromfile(filename::AbstractString, key::String) = fromfile(Units.defaultFloat, filename, key)