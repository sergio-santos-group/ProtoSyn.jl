module Sugars

using YAML

using ..ProtoSyn
using ..ProtoSyn.Builder


# resource directory for this module
const resource_dir = let
    modname = string(nameof(@__MODULE__))
    joinpath(ProtoSyn.resource_dir, modname)
end


"""
    grammar([T=Float64,] polyname::String)

Build a [`LGrammar`](@ref) for polysaccharide `polyname` from the `grammars.yml`
file available in the `Sugars` resource directory.
The returned L-grammar can then be used by the [`Builder.build`](@ref) function
to build the polymer.

# Examples
```julia-repl
julia> g = Sugars.grammar("amylose");
julia> pose = Builder.build(g, "AAA")
...
```
"""
function grammar(::Type{T}, key::String) where {T <: AbstractFloat}
    filename = joinpath(resource_dir, "grammars.yml")
    open(filename) do io
        @info "loading grammar for '$key' from file $filename"
        yml = YAML.load(io)
        lgfactory(T, yml[key])
    end
end

grammar(key::String) = grammar(Float64, key)

end
