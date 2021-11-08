"""
# TODO
"""
function scalling_exposed_only(Ωi::Union{Int, T}; hydrophobicity_map_value::T = 0.0, Ω::Union{Int, T} = 0.0) where {T <: AbstractFloat}
    if Ωi < Ω
        return T(0)
    else
        return hydrophobicity_map_value * (Ωi - Ω)
    end
end


"""
# TODO
"""
function non_scalling_exposed_only(Ωi::Union{Int, T}; hydrophobicity_map_value::T = 0.0, Ω::Union{Int, T} = 0.0) where {T <: AbstractFloat}
    if Ωi < Ω
        return T(0)
    else
        return hydrophobicity_map_value
    end
end


"""
# TODO
"""
function scalling_all_contributions(Ωi::Union{Int, T}; hydrophobicity_map_value::T = 0.0, Ω::Union{Int, T} = 0.0) where {T <: AbstractFloat}
    return hydrophobicity_map_value * (Ωi - Ω)
end


"""
# TODO
"""
function non_scalling_all_contributions(Ωi::Union{Int, T}; hydrophobicity_map_value::T = 0.0, Ω::Union{Int, T} = 0.0) where {T <: AbstractFloat}
    if Ωi < Ω
        return - hydrophobicity_map_value
    else
        return hydrophobicity_map_value
    end
end