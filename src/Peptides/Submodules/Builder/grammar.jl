export grammar

"""
    grammar([::Type{T}];[verbose::Bool = true]) where {T <: AbstractFloat}

Build a [`LGrammar`](@ref ProtoSyn.LGrammar) for peptides, taking as variables
the [`Fragment`](@ref ProtoSyn.Fragment) instances in the default resource
directory. If the optional type `T` is not provided, the `ProtoSyn.defaultFloat`
value will be used. The returned [`LGrammar`](@ref ProtoSyn.LGrammar) is
required for building peptides from [`Fragment`](@ref ProtoSyn.Fragment)
instances, for example. If `verbose` is set to `true` (is, by default), print
the loading status.

# Examples
```
julia> res_lib = ProtoSyn.Peptides.grammar()
```
"""
function grammar(::Type{T}) where {T <: AbstractFloat}
    filename = joinpath(Peptides.resource_dir, "grammars.yml")
    ProtoSyn.load_grammar_from_file(T, filename, "peptide")
end

grammar(;verbose::Bool = true) = grammar(ProtoSyn.Units.defaultFloat)
