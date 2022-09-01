export StochasticRule, LGrammar, Operation, Tautomer
export derive, build, lgfactory, load_grammar_from_file
export @seq_str

"""
    StochasticRule(p::T, rule::Pair{K, V}) where {T <: AbstractFloat, K, V}

Return a new [`StochasticRule`](@ref) instance with the given probability of
occurrence `p`. The `rule` is a `Pair{K, V}` where in most cases `K` is an
instance of type `String` (i.e.: a Key) and `V` is an instance of type
`Vector{String}` (i.e.: a Vector of instructions). These are also called of 
"production instructions" and define the result of deriving the given "key"
in any derivation. As an example, the pair `"A" => ["A", "ɑ", "A"]` would be
interpreted upon derivation, and the entry `A` would be _expanded_ to `AA`,
where both new `A` instances are joined by the `ɑ` operator, with a probability
of occurrence of `p`.

# Fields
* `p::T` - The probability of occurrence;
* `source::K` - The key of the `rule` on this [`StochasticRule`](@ref) instance;
* `production::V` - The resulting vector of the derivation of this [`StochasticRule`](@ref) instance on the given `source`.

# See also
[`LGrammar`](@ref)

# Examples
```jldoctest
julia> sr = StochasticRule(1.0, "A" => ["A", "ɑ", "A"])
A(p=1.0) -> ["A", "ɑ", "A"]
```
"""
struct StochasticRule{T <: AbstractFloat, K, V}
    p::T
    source::K
    production::V

    StochasticRule(p::T, rule::Pair{K, V}) where {T <: AbstractFloat, K, V} = begin
        (p > 1 || p < 0) && error("Invalid probability")
        new{T, K, V}(p, rule.first, rule.second)
    end
end

Base.show(io::IO, sr::StochasticRule) = begin
    println(io, sr.source, "(p=", sr.p, ") -> ", sr.production)
end

struct Operation
    op::Function
    r1::String
    r2::String
end

(op::Operation)(r1::Residue, f2::Union{Fragment, Pose}; residue_index = 1, segment_index = 1) = begin
    op.op(r1, f2, residue_index = residue_index, segment_index = segment_index)
end

function as_string(op::Operation)
    return "$(op.op) (Between atoms $(op.r1) & $(op.r2))"
end

Base.show(io::IO, op::Operation) = begin
    print(io, as_string(op))
end

struct Tautomer
    list::Vector{Fragment}
end

Base.getindex(t::ProtoSyn.Tautomer, i::Int64) = t.list[i]

function as_string(t::Tautomer)
    return "$(t.list[1])(And $(length(t.list)-1) other tautomer(s) available.)"
end

Base.show(io::IO, t::Tautomer) = begin
    print(io, as_string(t))
end

"""
    LGrammar{T <: AbstractFloat, K, V}(rules::Dict{K, Vector{StochasticRule{K,V}}}, variables::Dict{K, Fragment}, operators::Dict{K, Function}, defop::Opt{Function})

An [`LGrammar`](@ref) instance. Holds information regarding a stochastic
L-Grammar system, made up of a set of `variables` connectable by one or more
`operators`. Optionally, stochastic `rules` can randomly pick the operator to
apply, based on a set of weights.

    LGrammar{T, K, V}() where {T <: AbstractFloat, K, V}

Return an empty [`LGrammar`](@ref) instance.

# Fields:
* `rules::Dict{K, Vector{StochasticRule{K,V}}}` - A dictionary of [`StochasticRule`](@ref) instances indexed by the `variable` key over which the given rule will operate; 
* `variables::Dict{K, Fragment}` - A dictionary of variables ([`Fragment`](@ref) templates) indexed by the corresponding code;
* `operators::Dict{K, Function}` - A dictionary of operator `Function` instances indexed by a named `String`;
* `defop::Opt{Function}` - Default operator. If no operator is described in the given `derivation` (during the [`build`](@ref) process), uses this operator.

# See also
[`StochasticRule`](@ref) [`build`](@ref) [`load_grammar_from_file`](@ref)

!!! ukw "Note:"
    As a general rule, [`LGrammar`](@ref) instances are loaded from an .YML file
    (using the [`load_grammar_from_file`](@ref) method). Check this entry for a
    more in-depth description of the file format.

# Examples
```
julia> grammar = LGrammar{Float64, String, Vector{String}}()
LGrammar{Float64, String, Vector{String}}:
 Rules: None.
 Variables: None.
 Operators: None.

julia> grammar = ProtoSyn.Peptides.grammar
```
"""
mutable struct LGrammar{T <: AbstractFloat, K, V}
    rules::Dict{K, Vector{StochasticRule{T, K,V}}}
    variables::Dict{K, Union{Fragment, Tautomer}}
    operators::Dict{K, Operation}
    defop::Opt{Operation}
end

LGrammar{T, K, V}() where {T <: AbstractFloat, K, V} = LGrammar{T, K, V}(
    Dict{K, Vector{StochasticRule{T, K, V}}}(),
    Dict{K, Fragment}(),
    Dict{K, Operation}(),
    nothing
)

LGrammar() = LGrammar{Float64, String, Vector{String}}()

Base.show(io::IO, lgrammar::LGrammar{T, K, V}) where {T <: AbstractFloat, K, V} = begin
    print(io, "LGrammar{$T, $K, $V}:\n Rules:")
    if length(lgrammar.rules) == 0
        print(" None.")
    else
        for (key, rule) in lgrammar.rules
            print(io, "\n $key => $rule")
        end
    end

    print(io, "\n Variables:")
    if length(lgrammar.variables) == 0
        print(" None.")
    else
        for (key, variable) in lgrammar.variables
            print(io, "\n $key => $variable")
        end
    end

    print(io, "\n Operators:")
    if length(lgrammar.operators) == 0
        print(" None.")
    else
        for (key, operator) in lgrammar.operators
            print(io, "\n $key => $(as_string(operator))")
        end
    end
end


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
    !haskey(g.variables, k) && begin
        println(" - Key $k was not found in the provided LGrammar.")
        throw(KeyError(k))
    end

    if isa(g.variables[k], Tautomer)
        return g.variables[k][1]
    else
        return g.variables[k]
    end
end


getop(g::LGrammar, k) = begin
    !haskey(g.operators, k) && throw(KeyError(k))
    g.operators[k]
end

hasrule(g::LGrammar, r) = haskey(g.rules, r)

Base.push!(g::LGrammar{T, K, V}, r::StochasticRule{T, K,V}) where {T <: AbstractFloat, K, V} = begin
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

Base.setindex!(g::LGrammar{T, K, V}, x::Tautomer, key::K) where {T <: AbstractFloat, K, V} = begin
    g.variables[key] = x
    g
end

Base.setindex!(g::LGrammar{T, K, V}, x::Operation, key::K) where {T <: AbstractFloat, K, V} = begin
    g.operators[key] = x
    g
end

# ---
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
"""
    fragment(grammar::LGrammar{T, K, V}, derivation) where {T <: AbstractFloat, K, V}
    
Create and return a new [`Fragment`](@ref) ([`Pose`](@ref) instance with just a
single [`Segment`](@ref)) using the given `derivation` sequence on the provided
[`LGrammar`](@ref) `grammar` instructions. The main purpose of fragments is to
be temporary carriers of information, such as during the building process of a
new peptide from a sequence. Therefore, these structures often don't have any
real meaning and are, as such, deprived of a root/origin for the graph. Actual
structures should instead be of the slightly more complete type [`Pose`](@ref).

# See also
[`build`](@ref)

# Examples
```
julia> frag = fragment(res_lib, seq"AAA")
Fragment(Segment{/UNK:63875}, State{Float64}:
 Size: 30
 i2c: false | c2i: false
 Energy: Dict(:Total => Inf)
)
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
            push!(opstack, op.op)
        elseif letter == "["
            push!(stack, parent)
        elseif letter == "]"
            parent = pop!(stack)
        elseif isvar(grammar, letter)

            # ! This is a quick fix and should be changed in future versions
            if letter === "P"
                push!(opstack, grammar.operators["γ"].op)
            end

            frag = getvar(grammar, letter)

            frag2 = copy(frag)

            push!(seg, frag2.graph.items...) # Appending the residues to the segment
            append!(state, frag2.state)      # Merging the 2 states
            if parent !== nothing
                join = isempty(opstack) ? grammar.defop : pop!(opstack)
                join(parent.graph[end], frag2) # Adding ascendents and bonds correctly
            end
            parent = frag2
        else
            @warn "Residue $letter in fragment derivation can't be identified either as a variable or as an operation."
        end
    end

    reindex(seg)
    seg.id = state.id = ProtoSyn.genid()

    return Pose(seg, state)
end

#endregion fragment


"""
    opfactory(args::Any)

Return the `operation` function (as a closure) given the input arguments `args`
(normally read from a grammar file). 

# See also
[`lgfactory`](@ref)

!!! ukw "Note:"
    The resulting operation function is responsible for setting the internal
    coordinates of residues in the system when connecting, building and
    manipulating poses.

!!! ukw "Note:"
    This is an internal method of ProtoSyn and shouldn't normally be used
    directly.
"""
function opfactory(args::Any)
    # Note: residue_index is the index on the fragment/pose ('f2') used to
    # actually connect to 'r1'
    op = function(r1::Residue, f2::Union{Fragment, Pose}; residue_index = 1, segment_index = 1)
        # residue_index is of f2
        if typeof(f2) == Fragment
            r2 = f2.graph[residue_index]
        else
            r2 = f2.graph[segment_index][residue_index]
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

        if haskey(args, "dependents")
            for (atname, entry) in args["dependents"]
                !(atname in r2) && continue # Proline has no H
                atomstate = state[r2[atname]]
                parent = state[r2[atname].parent]

                for (key, dependence) in entry
                    !haskey(dependence, "measure") && begin
                        error("No 'measure' key found in dependence record")
                    end
                    measured_atom = dependence["measure"]
                    r1[measured_atom] === nothing && continue
                    value = state[r1[measured_atom]].ϕ + state[r1[measured_atom].parent].Δϕ

                    haskey(dependence, "offset") && begin
                        value += dependence["offset"]
                    end

                    setproperty!(atomstate, Symbol(key), value - parent.Δϕ)
                end
            end
        end
    end

    return Operation(op, args["residue1"], args["residue2"])
end



"""
    lgfactory([::Type{T}], template::Dict) where {T <: AbstractFloat}

Create an `LGrammar` instance from the contencts of a `template` Dict (normally
read from a grammar file). Any numerical entry is parsed to the provided type
`T` (or `Units.defaultFloat` if no type is provided). The `operators` entry is
parsed by the [`opfactory`](@ref) method. Return the parsed [`LGrammar`](@ref)
instance.

# See also
[`LGrammar`](@ref) [`load_grammar_from_file`](@ref) [`opfactory`](@ref)

!!! ukw "Note:"
    This is an internal method of ProtoSyn and shouldn't normally be used
    directly.
"""
function lgfactory(::Type{T}, template::Dict) where T
    grammar = LGrammar{T, String, Vector{String}}()

    vars = template["variables"]
    for (key, name) in vars
        if isa(name, Vector{String})
            poses = Vector{Fragment}()
            for (i, _name) in enumerate(name)
                filename = joinpath(ProtoSyn.resource_dir, _name)
                @info "Loading variable '$key' (tautomer $i) from $filename"
                pose = ProtoSyn.load(T, filename)
                pose.graph.name = pose.graph[1].name # ! Hack for long filenames
                push!(poses, fragment(pose))
            end
            grammar[key] = Tautomer(poses)
        else
            filename = joinpath(ProtoSyn.resource_dir, name)
            @info "Loading variable '$key' from $filename"
            pose = ProtoSyn.load(T, filename)
            pose.graph.name = pose.graph[1].name # ! Hack for long filenames
            grammar[key] = fragment(pose)
        end

    end

    if haskey(template, "operators")
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

            if haskey(opargs, "dependents")
                dependents = opargs["dependents"]
                for (k, v) in dependents
                    for (_k, _v) in v
                        !haskey(dependents[k][_k], "offset") && continue
                        dependents[k][_k]["offset"] = ProtoSyn.Units.tonumber(T, _v["offset"])
                    end
                end
            end

            @info "Loading operator $opname"
            grammar[opname] = opfactory(opargs)
        end
    end
    
    if haskey(template, "defop")
        grammar.defop = getop(grammar, template["defop"])
    end

    if haskey(template, "rules")
        for (key, rules) in template["rules"]
            @info "Loading productions for rule $key"
            for rule in rules
                sr = StochasticRule(rule["p"], key => rule["production"])
                push!(grammar, sr)
                @info "  $sr"
            end
        end
    end

    return grammar
end

lgfactory(template::Dict) = lgfactory(Float64, template)


"""
    load_grammar_from_file([::Type{T}], filename::AbstractString, key::String) where {T <: AbstractFloat}

Create an [`LGrammar`](@ref) instance from the contents of a grammar file (in
.YML format) under the `key` entry. The file contents are parsed by the
[`lgfactory`](@ref) method. Any numerical entry is parsed to the provided type
`T` (or `Units.defaultFloat` if no type is provided). Return the parsed
[`LGrammar`](@ref) instance. automatically calls
[`load_grammar_extras_from_file!`](@ref).

# See also
[`LGrammar`](@ref) [`lgfactory`](@ref)

# Examples
```
julia> lgrammar = load_grammar_from_file(Float64, filename, "peptide")

julia> lgrammar = load_grammar_from_file(filename, "peptide")
```
"""
function load_grammar_from_file(::Type{T}, filename::AbstractString, key::String) where {T <: AbstractFloat}
    ProtoSyn.load_grammar_extras_from_file!(T, filename, key)
    open(filename) do io
        @info "loading grammar from file $filename"
        yml = YAML.load(io)
        return lgfactory(T, yml[key])
    end
end

load_grammar_from_file(filename::AbstractString, key::String) = begin
    return load_grammar_from_file(Units.defaultFloat, filename, key)
end


"""
    load_grammar_extras_from_file!([::Type{T}], filename::AbstractString, key::String) where {T <: AbstractFloat}

Loads the `key` entry in the given [`LGrammar`](@ref) .YML file (`filename`)
extras into the correct global variables in ProtoSyn. Any numerical entry is
parsed to the provided type `T` (or `Units.defaultFloat` if no type is
provided).

The extra info loaded by this method is:
+ Any `alt` entry is added to `ProtoSyn.alt_residue_names`

!!! ukw "Note:"
    Other modules (such as [Peptides](@ref)) may retrieve extra information from the [`LGrammar`](@ref) file. As such, these modules often include an expanded method for [`load_grammar_extras_from_file!`](@ref).

!!! ukw "Note:"
    This method is automatically called from [`load_grammar_from_file`](@ref). This is the recommended way to load an [`LGrammar`](@ref) (this method shouldn't be called as a standalone for most applications).

# Examples
```
julia> ProtoSyn.load_grammar_extras_from_file!(ProtoSyn.resource_dir*"/Peptides/grammars.yml", "default")
```
"""
function load_grammar_extras_from_file!(::Type{T}, filename::AbstractString, key::String) where {T <: AbstractFloat}
    
    function load_alt_names!(var_yml::Dict{Any, Any})
        if !(var_yml["name"] in keys(ProtoSyn.alt_residue_names))
            ProtoSyn.alt_residue_names[var_yml["name"]] = Vector{String}([var_yml["name"]])
        end
        if "alt" in keys(var_yml)
            for alt in var_yml["alt"]
                if !(alt in ProtoSyn.alt_residue_names[var_yml["name"]])
                    push!(ProtoSyn.alt_residue_names[var_yml["name"]], alt)
                end

                if alt in keys(ProtoSyn.alt_residue_names)
                    if !(var_yml["name"] in ProtoSyn.alt_residue_names[alt])
                        push!(ProtoSyn.alt_residue_names[alt], var_yml["name"])
                    end
                else
                    ProtoSyn.alt_residue_names[alt] = Vector{String}([alt, var_yml["name"]])
                end
            end
        end
    end

    # Re-read the YML grammar file
    yml = ProtoSyn.read_yml(filename)[key]

    vars = yml["variables"]
    for (key, name) in vars
        if isa(name, Vector{String})
            for _name in name
                var_filename = joinpath(ProtoSyn.resource_dir, _name)
                var_yml = ProtoSyn.read_yml(var_filename)
                load_alt_names!(var_yml)
            end
        else
            var_filename = joinpath(ProtoSyn.resource_dir, name)
            var_yml = ProtoSyn.read_yml(var_filename)
            load_alt_names!(var_yml)
        end
    end
end

load_grammar_extras_from_file!(filename::AbstractString, key::String) = begin
    ProtoSyn.load_grammar_extras_from_file!(ProtoSyn.Units.defaultFloat, filename, key)
end

function expand(grammar::LGrammar{T}, derivation; n::Int = 1) where {T <: AbstractFloat}
    for j in 1:n
        new_derivation = []
        for i in derivation
            new_derivation = vcat(new_derivation, ProtoSyn.getrule(grammar, i))
        end

        derivation = copy(new_derivation)
    end

    return derivation
end