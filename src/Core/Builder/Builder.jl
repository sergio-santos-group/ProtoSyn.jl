module Builder

using YAML

using ..ProtoSyn
using ..ProtoSyn.Units: tonumber

export StochasticRule, LGrammar
export derive, build, lgfactory, fromfile
export @seq_str

struct StochasticRule{K,V}
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

mutable struct LGrammar{K,V}
    rules::Dict{K, Vector{StochasticRule{K,V}}}
    variables::Dict{K, Fragment}
    operators::Dict{K, Function}
    defop::Opt{Function}
end


LGrammar{K,V}() where {K,V} = LGrammar(
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


function fragment(::Type{T}, grammar::LGrammar, derivation) where {T<:AbstractFloat}

    state = State{T}()
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
fragment(grammar::LGrammar, derivation) = fragment(Float64, grammar, derivation)


"""
    build([T=Float64,] grammar::LGrammar, derivation)

Build a `Pose{Topology}`.
"""
function build(::Type{T}, grammar::LGrammar, derivation) where {T<:AbstractFloat}
    top = Topology("UNK", 1)
    state = State{T}()
    state.id = top.id
    pose = Pose(top, state)

    if !isempty(derivation)
        frag = fragment(T, grammar, derivation)
        append!(pose, frag) # Appending the fragment (which is a segment) to the Topology
        
        ProtoSyn.request_i2c(state; all=true)
    end
    pose
end
build(grammar::LGrammar, derivation) = build(Float64, grammar, derivation)


"""
    @seq_str -> Vector{String}

Construct a vector of strings from the provided string. 

# Examples
```jldoctest
julia> seq"ABC"
3-element Array{String,1}:
 "A"
 "B"
 "C"
```
"""
macro seq_str(s); [string(c) for c in s]; end

"""
    opfactory(args) -> Function

Closure
"""
function opfactory(args)
    # Note: residue_index is the index on the fragment/pose ('f2') used to
    # actually connect to 'r1'
    return function(r1::Residue, f2::Union{Fragment, Pose}; residue_index = 1)
        # Index is of f2
        if typeof(f2) == Fragment
            r2 = f2.graph[residue_index]
        else
            r2 = f2.graph[1][residue_index]
        end
        ProtoSyn.join(r1, args["residue1"], r2, args["residue2"]) # Connects specifically C to N (in case of proteins) -> Adds bonds and sets parents
        state = f2.state
        
        if haskey(args, "presets")
            for (atname, presets) in args["presets"]
                # println(" Atom $(r2[atname])")
                atomstate = state[r2[atname]]
                for (key, value) in presets
                    setproperty!(atomstate, Symbol(key), value)
                end
            end
        end
        return
        
        if haskey(args, "offsets")
            for (atname, offset) in args["offsets"]
                setoffset!(state, r2[atname], offset)
            end
        end
    end
end



"""
    lgfactory([T=Float64,] template::Dict) -> LGrammar{String,Vector{String}}


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
    grammar = LGrammar{String,Vector{String}}()

    vars = template["variables"]
    for (key, name) in vars
        filename = joinpath(ProtoSyn.resource_dir, name)
        @info "Loading variable '$key' from $filename"
        pose = ProtoSyn.load(T, filename)
        grammar[key] = ProtoSyn.fragment(pose)
    end

    ops = template["operators"]
    for (opname, opargs) in ops
        if haskey(opargs, "presets")
            for presets in values(opargs["presets"])
                for (k, v) in presets
                    presets[k] = tonumber(v)
                end
            end
        end

        if haskey(opargs, "offsets")
            offsets = opargs["offsets"]
            for (k,v) in offsets
                offsets[k] = tonumber(v)
            end
        end

        @info "Loading operator $opname"
        grammar[opname] = opfactory(opargs)
    end
    grammar.defop = Builder.getop(grammar, template["defop"])

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


@inline function bond(pose::Pose{Topology}, at1::Atom, at2::Atom, grammar::Builder.LGrammar; op = "α")

    @assert ProtoSyn.segment(at1) === ProtoSyn.segment(at2) "can only bond atoms within the same segment"

    if at2.parent == ProtoSyn.origin(pose.graph)
        popparent!(at2)
        popparent!(at2.container)
    end

    grammar.operators[op](at1.container, pose, residue_index = at2.container.index)
    reindex(pose.graph)
    ProtoSyn.request_i2c(pose.state)

    return pose
end


function fromfile(::Type{T}, filename::AbstractString, key::String) where T
    open(filename) do io
        @info "loading grammar from file $filename"
        yml = YAML.load(io)
        lgfactory(T, yml[key])
    end
end


function append_residues!(pose::Pose{Topology}, residue::Residue, grammar::LGrammar, derivation; op = "α")

    frag = Builder.fragment(grammar, derivation)

    # Inserts the fragment residues in the Pose graph and sets residue.container
    # println("Insert residues in position $(residue.index + 1)/$(length(residue.container.items))")
    insert!(residue.container, residue.index + 1, frag.graph.items)

    # Inserts the fragment atoms state in the pose.state
    # println("Insert residues in position $(residue.items[end].index + 1)/$(length(pose.state.items)-3)")
    insert!(pose.state, residue.items[end].index + 1, frag.state)

    # Sets:
    # - Distance/angle/dihedrals in the fragment first residue
    # - Parent/children in the newly bonded atoms
    # - Parent/children in the newly bonded residues
    # - Atom bonds in the newly bonded atoms

    reindex(pose.graph, set_ascendents = false)

    grammar.operators[op](residue, pose, residue_index = residue.index + 1)

    # Reindex and define new ascendents
    reindex(pose.graph)
    ProtoSyn.request_i2c(pose.state; all=true)
end


function insert_residues!(pose::Pose{Topology}, residue::Residue, grammar::LGrammar, derivation; op = "α", connect_upstream = true)

    residue_index = residue.index
    frag = Builder.fragment(grammar, derivation)

    insert!(residue.container, residue.index, frag.graph.items)
    insert!(pose.state, residue.items[1].index, frag.state)

    # Remove all references of the origin as the parent of any of its children
    # and remove said children from origin.children
    pose_origin = ProtoSyn.origin(pose.graph) # This is an Atom
    parent_is_origin = residue.parent == pose_origin.container
    if parent_is_origin
        for child in pose_origin.children
            child in residue.items && popparent!(child)
        end
        
        # Add the first atom of the appendage as a child of origin AND set its
        # parent as being the origin
        _root = ProtoSyn.root(frag.graph)
        setparent!(_root, pose_origin)
        setparent!(_root.container, pose_origin.container)
        popparent!(residue)
    end

    # Reindex to account for the inserted residues
    reindex(pose.graph, set_ascendents = false)
    grammar.operators[op](frag.graph[end], pose, residue_index = residue_index + length(frag.graph))


    # Case we are inserting between two pre-existing residues (so far, the same
    # operation will be used in both cases. Might be useful to differentiate
    # between left and right operation.)
    if connect_upstream
        for child in pose_origin.children
            child in residue.container[residue_index].items && popparent!(child)
        end
        popparent!(residue.container[residue_index])

        grammar.operators[op](residue.container[residue_index - 1], pose, residue_index = residue_index)
    end
    
    # Reindex to set correct ascendents
    reindex(pose.graph)
    ProtoSyn.request_i2c(pose.state; all=true)
end


# function mutate!(pose::Pose{Topology}, residue::Residue, grammar::LGrammar, derivation; op = "α")

#     @assert length(derivation) == 1 "Derivation must have length = 1."
#     residue_index = residue.index

#     insert_residues!(pose, residue, grammar, derivation; op = op)
#     # pop!(pose, residue.container[residue.index])



#     # bond(pose, residue.container[residue_index], residue.container[residue_index + 1], grammar, op = op)
# end

end