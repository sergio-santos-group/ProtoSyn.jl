using YAML

export grammar

"""
# TODO
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
function load_grammar_from_file(::Type{T}, filename::AbstractString, key::String) where {T <: AbstractFloat}
    filename = joinpath(Peptides.resource_dir, filename)
    grammar = ProtoSyn.load_grammar_from_file(T, filename, key)
    load_grammar_extras_from_file!(T, filename, key)

    return grammar
end

load_grammar_from_file(filename::AbstractString, key::String) = begin
    load_grammar_from_file(ProtoSyn.Units.defaultFloat, filename, key)
end

"""
# TODO
"""
function load_grammar_extras_from_file!(::Type{T}, filename::AbstractString, key::String) where {T <: AbstractFloat}
    
    function read_yml(_filename::String)
        open(_filename, "r") do io
            return YAML.load(io)
        end
    end

    # Re-read the YML grammar file
    yml = read_yml(filename)[key]
    
    vars = yml["variables"]
    for (key, name) in vars
        var_filename = joinpath(ProtoSyn.resource_dir, name)
        var_yml = read_yml(var_filename)
        
        # 1. Add chi entries to Peptides.Dihedral.chi_dict
        if "chis" in keys(var_yml)
            Peptides.Dihedral.chi_dict[var_yml["name"]] = var_yml["chis"]
        end

        # 2. Add entries to the Peptides.available_aminoacids
        Peptides.available_aminoacids[var_yml["code"][1]] = true

        # 3. Add entries to Peptides.three_2_one & Peptides.one_2_three
        Peptides.three_2_one[var_yml["name"]] = var_yml["code"][1]
        if "alt" in keys(var_yml)
            for alt_name in var_yml["alt"]
                Peptides.three_2_one[alt_name] = var_yml["code"][1]
            end
        end
        Peptides.one_2_three[var_yml["code"][1]] = var_yml["name"]
    end
end

load_grammar_extras_from_file!(filename::AbstractString, key::String) = begin
    load_grammar_extras_from_file!(ProtoSyn.Units.defaultFloat, filename, key)
end


"""
# TODO
"""
load_default_grammar(::Type{T}) where {T <: AbstractFloat} = load_grammar_from_file(T, "grammars.yml", "default")
load_default_grammar() = load_default_grammar(ProtoSyn.Units.defaultFloat)


"""
# TODO
"""
grammar = load_default_grammar()