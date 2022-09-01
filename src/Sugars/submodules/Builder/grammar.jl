export grammar

using ProtoSyn

# Resource directory for this module
const resource_dir = let
    modname = string(nameof(@__MODULE__))
    joinpath(ProtoSyn.resource_dir, modname)
end


"""
    load_default_grammar([::Type{T}, polyname::String]) where {T <: AbstractFloat}

Build a [`LGrammar`](@ref ProtoSyn.LGrammar) for polysaccharide `polyname` from
the `grammars.yml` file available in the `Sugars` resource directory ("amylose",
by default). The returned [`LGrammar`](@ref ProtoSyn.LGrammar) can then be used
by the [`ProtoSyn.build`](@ref) function to build the polymer.

# Examples
```
julia> g = Sugars.grammar;

julia> pose = ProtoSyn.build(g, seq"AAAβB[ɣCɣCɣC]AAA")
```
"""
function load_default_grammar(::Type{T}, polyname::String) where {T <: AbstractFloat}
    filename = joinpath(resource_dir, "grammars.yml")
    ProtoSyn.load_grammar_from_file(T, filename, polyname)
end

grammar = load_default_grammar(ProtoSyn.Units.defaultFloat, "amylose")