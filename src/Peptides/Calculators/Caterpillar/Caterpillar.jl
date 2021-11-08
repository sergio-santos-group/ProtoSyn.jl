module Caterpillar

    using ProtoSyn
    using ProtoSyn.Calculators: EnergyFunctionComponent, full_distance_matrix
    using LinearAlgebra


    # Burial degree calculation
    include("neighbour_count.jl")
    include("neighbour_vector.jl")

    # Identification curve
    include("identification_curve.jl")

    # Hydrophobicity Weight
    include("hydrophobicity_weight.jl")

    function get_default_solvation_energy(;α::T = 1.0) where {T <: AbstractFloat}
        return EnergyFunctionComponent(
            "Caterpillar_Solvation",
            neighbour_vector,
            Dict{Symbol, Any}(
                :identification_curve  => sigmoid,
                :selection             => an"CB" | (an"HA2" & rn"GLY"),
                :rmax                  => 9.0,
                :sc                    => 1.0,
                :Ω                     => 24.0,
                :hydrophob_map         => ProtoSyn.Peptides.doolitle_hydrophobicity_mod7,
                :hydrophobicity_weight => scalling_all_contributions),
            α,
            false)
    end

end