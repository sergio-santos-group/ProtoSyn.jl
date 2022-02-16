module Caterpillar

    using ProtoSyn
    using ProtoSyn.Calculators: EnergyFunctionComponent, full_distance_matrix
    using LinearAlgebra


    # Burial degree calculation
    null_identification_curve(x::Float64; rmax::Float64 = 0.0, sc::Float64 = 0.0) = 0.0
    null_hydrophobicity_weight(x::Float64; hydrophobicity_map_value::Float64 = 0.0, Ω::Float64 = 0.0) = 0.0

    include("neighbour_count.jl")
    include("neighbour_vector.jl")

    # Identification curve
    include("identification_curve.jl")

    # Hydrophobicity Weight
    include("hydrophobicity_weight.jl")

    function get_default_solvation_energy(;α::T = 1.0) where {T <: AbstractFloat}
        return EnergyFunctionComponent(
            "Caterpillar_Solv",
            neighbour_vector,
            an"CB" | (an"HA2" & rn"GLY"),
            Dict{Symbol, Any}(
                :identification_curve  => sigmoid,
                :rmax                  => 9.0,
                :sc                    => 1.0,
                :Ω                     => 24.0,
                :hydrophobicity_map    => ProtoSyn.Peptides.doolitle_hydrophobicity_mod7,
                :hydrophobicity_weight => nv_scalling_all_contributions),
            α,
            false)
    end

end