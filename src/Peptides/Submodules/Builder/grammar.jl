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
function join_grammars!(g1::LGrammar{T, K, V}, g2::LGrammar{T, K, V}) where {T <: AbstractFloat, K, V}
    for (key, value) in g2.rules
        if key in keys(g1.rules)
            @warn "Rule $key found in existing grammar with same key. Skipping addition of rule from new grammar. Check if this is the desired behaviour."
        else
            g1.rules[key] = value
        end
    end

    for (key, value) in g2.variables
        if key in keys(g1.variables)
            @warn "Variable $key found in existing grammar with same key. Skipping addition of variable from new grammar. Check if this is the desired behaviour."
        else
            g1.variables[key] = value
        end
    end

    for (key, value) in g2.operators
        if key in keys(g1.operators)
            @warn "Operator $key found in existing grammar with same key. Skipping addition of operator from new grammar. Check if this is the desired behaviour."
        else
            g1.operators[key] = value
        end
    end

    @warn "Default operator of existing grammar was not changed. Skipping addition of default operator from new grammar. Check if this is the desired behaviour."

    return g1
end

"""
# TODO
"""
function load_grammar_from_file!(::Type{T}, grammar::LGrammar{T, K, V}, filename::AbstractString, key::String) where {T <: AbstractFloat, K, V}
    new_grammar = load_grammar_from_file(T, filename, key)
    join_grammars!(grammar, new_grammar)

    return grammar
end

load_grammar_from_file!(grammar::LGrammar{T, K, V}, filename::AbstractString, key::String) where {T <: AbstractFloat, K, V} = begin
    load_grammar_from_file!(ProtoSyn.Units.defaultFloat, grammar, filename, key)
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
        Peptides.one_2_three[var_yml["code"][1]] = var_yml["name"]

        # 4. Add same entries for any alternative name
        if "alt" in keys(var_yml)
            for alt_name in var_yml["alt"]
                Peptides.three_2_one[alt_name] = var_yml["code"][1]
                Peptides.Dihedral.chi_dict[alt_name] = var_yml["chis"]
            end
        end
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