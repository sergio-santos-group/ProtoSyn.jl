using YAML

export grammar

"""
    load_grammar_from_file([::Type{T}], filename::AbstractString, key::String) where {T <: AbstractFloat}

Load a [`LGrammar`](@ref ProtoSyn.LGrammar) from the given `filename` (this
should in .yml format, only loads the given `key`), typed to be of type
`T <: AbstractFloat`. If not provided, will use `ProtoSyn.Units.defaultFloat`.

# See also
[`load_grammar_from_file!`](@ref)

# Examples
```
julia> ProtoSyn.Peptides.load_grammar_from_file("grammars.yml", "default")
LGrammar{Float64, String, Vector{String}}:
 (...)
```
"""
function load_grammar_from_file(::Type{T}, filename::AbstractString, key::String) where {T <: AbstractFloat}
    filename = joinpath(Peptides.resource_dir, filename)
    grammar = ProtoSyn.load_grammar_from_file(T, filename, key)
    Peptides.load_grammar_extras_from_file!(T, filename, key)

    return grammar
end

load_grammar_from_file(filename::AbstractString, key::String) = begin
    load_grammar_from_file(ProtoSyn.Units.defaultFloat, filename, key)
end


"""
    join_grammars!(g1::LGrammar{T, K, V}, g2::LGrammar{T, K, V}) where {T <: AbstractFloat, K, V}

Add [`LGrammar`](@ref) `g2` rules, variables and operators to [`LGrammar`](@ref)
`g1`, joining both instances.

# Examples
```
julia> g1 = ProtoSyn.Peptides.load_grammar_from_file("grammars.yml", "default");

julia> g2 = ProtoSyn.Peptides.load_grammar_from_file("grammars.yml", "ncaa");

julia> ProtoSyn.Peptides.join_grammars!(g1, g2)
LGrammar{Float64, String, Vector{String}}:
 (...)
```
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
    load_grammar_from_file!([::Type{T}], grammar::LGrammar{T, K, V}, filename::AbstractString, key::String) where {T <: AbstractFloat, K, V}

Load a [`LGrammar`](@ref ProtoSyn.LGrammar) from the given `filename` (this
should in .yml format, only loads the given `key`), typed to be of type
`T <: AbstractFloat`. If not provided, will use `ProtoSyn.Units.defaultFloat`.
Automatically add the loaded [`LGrammar`](@ref) instance to the given `grammar`
(using the [`join_grammars!`](@ref) method).

# See also
[`load_grammar_from_file`](@ref)

# Examples
```
julia> g1 = ProtoSyn.Peptides.load_grammar_from_file("grammars.yml", "default");

julia> g2 = ProtoSyn.Peptides.load_grammar_from_file!(g1, "grammars.yml", "ncaa")
LGrammar{Float64, String, Vector{String}}:
 (...)
```
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
    load_grammar_extras_from_file!([::Type{T}], filename::AbstractString, key::String) where {T <: AbstractFloat}

Fill missing information in the following ProtoSyn constants with extra
details from the loaded [`LGrammar`](@ref) .yml `filename` (under the `key`
entry, typed as `T <: AbstractFloat`, if not provided will use
`ProtoSyn.Units.defaultFloat`):

* `ProtoSyn.Peptides.chi_dict`
* `ProtoSyn.Peptides.available_aminoacids`
* `ProtoSyn.three_2_one`
* `ProtoSyn.one_2_three`

# See also
[`load_grammar_from_file`](@ref)

# Examples
```
julia> ProtoSyn.Peptides.load_grammar_extras_from_file!("grammars.yml", "default")
```
"""
function load_grammar_extras_from_file!(::Type{T}, filename::AbstractString, key::String) where {T <: AbstractFloat}

    function load_extras!(var_yml::Dict{Any, Any})
        # 1. Add chi entries to Peptides.chi_dict
        if "chis" in keys(var_yml)
            Peptides.chi_dict[var_yml["name"]] = var_yml["chis"]
        end

        # 2. Add entries to the Peptides.available_aminoacids
        Peptides.available_aminoacids[var_yml["code"][1]] = true

        # 3. Add entries to ProtoSyn.three_2_one & ProtoSyn.one_2_three
        ProtoSyn.three_2_one[var_yml["name"]] = var_yml["code"][1]
        if var_yml["code"][1] in keys(ProtoSyn.one_2_three)
            @info "Variable $(var_yml["code"][1]) already in the Protosyn.three_2_one dictionary (Currently assigned to $(ProtoSyn.one_2_three[var_yml["code"][1]]))."
        else
            ProtoSyn.one_2_three[var_yml["code"][1]] = var_yml["name"]
        end

        # 4. Add same entries for any alternative name
        if "alt" in keys(var_yml)
            for alt_name in var_yml["alt"]
                ProtoSyn.three_2_one[alt_name] = var_yml["code"][1]
                if alt_name in keys(Peptides.chi_dict)
                    @info "Variable $(alt_name) already in the Peptides.chi_dict (Currently assigned to $(Peptides.chi_dict[alt_name])."
                else
                    Peptides.chi_dict[alt_name] = var_yml["chis"]
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
                load_extras!(var_yml)
            end
        else
            var_filename = joinpath(ProtoSyn.resource_dir, name)
            var_yml = ProtoSyn.read_yml(var_filename)
            load_extras!(var_yml)
        end
    end
end

load_grammar_extras_from_file!(filename::AbstractString, key::String) = begin
    Peptides.load_grammar_extras_from_file!(ProtoSyn.Units.defaultFloat, filename, key)
end


"""
    load_default_grammar([::Type{T}])

Load the default Peptides [`LGrammar`](@ref) (from the default
resources/Peptides/ directory, types as `T <: AbstractFloat`, if not provided
will use `ProtoSyn.Units.defaultFloat`).

!!! ukw "Note:"
    The default Peptides [`LGrammar`](@ref) is automatically loaded when using ProtoSyn. It can be found at `Peptides.grammar`.

# Examples
```
julia> ProtoSyn.Peptides.load_default_grammar()
LGrammar{Float64, String, Vector{String}}:
 (...)
```
"""
load_default_grammar(::Type{T}) where {T <: AbstractFloat} = load_grammar_from_file(T, "grammars.yml", "default")
load_default_grammar() = load_default_grammar(ProtoSyn.Units.defaultFloat)

@info " | Loading default peptides grammar"

grammar = load_default_grammar()