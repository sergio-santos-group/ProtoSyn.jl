module Sugars

using ..ProtoSyn
using ..ProtoSyn

export grammar

# resource directory for this module
const resource_dir = let
    modname = string(nameof(@__MODULE__))
    joinpath(ProtoSyn.resource_dir, modname)
end


"""
    grammar([T=Float64,] polyname::String)

Build a [`LGrammar`](@ref) for polysaccharide `polyname` from the `grammars.yml`
file available in the `Sugars` resource directory.
The returned L-grammar can then be used by the [`ProtoSyn.build`](@ref) function
to build the polymer.

# Examples
```julia-repl
julia> g = Sugars.grammar("amylose");
julia> pose = ProtoSyn.build(g, "AAA")
...
```
"""
function grammar(::Type{T}, polyname::String) where {T <: AbstractFloat}
    filename = joinpath(resource_dir, "grammars.yml")
    ProtoSyn.load_grammar_from_file(T, filename, polyname)
end
grammar(polyname::String) = grammar(Float64, polyname)

end
