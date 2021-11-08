"""
    scalling_exposed_only(Ωi::Union{Int, T}; hydrophobicity_map_value::T = 0.0, Ω::Union{Int, T} = 0.0) where {T <: AbstractFloat}

Return the hydrophobicity weight `w2` where all buried residues contributions
are set to 0.0. The exposed residues contribution is a scaled difference between
the given `Ωi` and the cut-off `Ω`, multiplied by the
`hydrophobicity_map_value`.

# See also
[`non_scalling_exposed_only`](@ref) [`scalling_all_contributions`](@ref)
[`non_scalling_all_contributions`](@ref)

# Examples
```
julia> ProtoSyn.Peptides.Calculators.Caterpillar.scalling_exposed_only(10.0, hydrophobicity_map_value = 1.0, Ω = 5.0)
5.0
```
"""
function scalling_exposed_only(Ωi::Union{Int, T}; hydrophobicity_map_value::T = 0.0, Ω::Union{Int, T} = 0.0) where {T <: AbstractFloat}
    if Ωi < Ω
        return T(0)
    else
        return hydrophobicity_map_value * (Ωi - Ω)
    end
end


"""
    non_scalling_exposed_only(Ωi::Union{Int, T}; hydrophobicity_map_value::T = 0.0, Ω::Union{Int, T} = 0.0) where {T <: AbstractFloat}

Return the hydrophobicity weight `w2` where all buried residues contributions
are set to 0.0. The exposed residues contribution is the non scaled
`hydrophobicity_map_value`.

# See also
[`scalling_exposed_only`](@ref) [`scalling_all_contributions`](@ref)
[`non_scalling_all_contributions`](@ref)

# Examples
```
julia> ProtoSyn.Peptides.Calculators.Caterpillar.non_scalling_exposed_only(10.0, hydrophobicity_map_value = 1.0, Ω = 5.0)
1.0
```
"""
function non_scalling_exposed_only(Ωi::Union{Int, T}; hydrophobicity_map_value::T = 0.0, Ω::Union{Int, T} = 0.0) where {T <: AbstractFloat}
    if Ωi < Ω
        return T(0)
    else
        return hydrophobicity_map_value
    end
end


"""
    scalling_all_contributions(Ωi::Union{Int, T}; hydrophobicity_map_value::T = 0.0, Ω::Union{Int, T} = 0.0) where {T <: AbstractFloat}

Return the hydrophobicity weight `w2` where all residue contributions are
considered. Both the exposed and buried residues contributions are a scaled
difference between the given `Ωi` and the cut-off `Ω`, multiplied by the
`hydrophobicity_map_value`.

# See also
[`scalling_exposed_only`](@ref) [`non_scalling_exposed_only`](@ref)
[`non_scalling_all_contributions`](@ref)

# Examples
```
julia> ProtoSyn.Peptides.Calculators.Caterpillar.scalling_all_contributions(10.0, hydrophobicity_map_value = 1.0, Ω = 20.0)
-10.0
```
"""
function scalling_all_contributions(Ωi::Union{Int, T}; hydrophobicity_map_value::T = 0.0, Ω::Union{Int, T} = 0.0) where {T <: AbstractFloat}
    return hydrophobicity_map_value * (Ωi - Ω)
end


"""
    non_scalling_all_contributions(Ωi::Union{Int, T}; hydrophobicity_map_value::T = 0.0, Ω::Union{Int, T} = 0.0) where {T <: AbstractFloat}

Return the hydrophobicity weight `w2` where all residue contributions are
considered. Both the exposed and buried residues contributions are the non
scaled `hydrophobicity_map_value`.

# See also
[`scalling_exposed_only`](@ref) [`non_scalling_exposed_only`](@ref)
[`scalling_all_contributions`](@ref)

# Examples
```
julia> ProtoSyn.Peptides.Calculators.Caterpillar.non_scalling_all_contributions(10.0, hydrophobicity_map_value = 1.0, Ω = 20.0)
-1.0
```
"""
function non_scalling_all_contributions(Ωi::Union{Int, T}; hydrophobicity_map_value::T = 0.0, Ω::Union{Int, T} = 0.0) where {T <: AbstractFloat}
    if Ωi < Ω
        return - hydrophobicity_map_value
    else
        return hydrophobicity_map_value
    end
end