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

    """
        get_default_caterpillar_solvation_energy(;[α::T = 1.0]) where {T <: AbstractFloat}

    Return the default Caterpillar solvation energy energy [`EnergyFunctionComponent`](@ref).
    `α` sets the component weight (on an
    [`EnergyFunction`](@ref ProtoSyn.Calculators.EnergyFunction) instance, `1.0`
    by default). This function employs [`neighbour_vector`](@ref) as the `:calc`
    function.
    
    # Settings
    * `hydrophobicity_weight::Function` - Returns the hydrophobicity weight of a given burial degree based on a hydrophobicity map value and default buried state;
    * `identification_curve::Function` - Returns the burial degree weight given an inter-atomic distance, rmax and slope control value;
    * `Ω::Union{Int, T}` - The default buried state (the average, buried degrees above this value are considered buried);
    * `sc::T = 1.0` - Slope control, used in sigmoid identification curves to control how sharp to consider an interaction as buried or not.
    * `hydrophobicity_map::Dict{String, T}` - A 1-on-1 map between [`Residue`](@ref) types and hydrophobicity indexes;
    * `rmax::T` - The cutoff value to consider some inter-atomic interaction as buriable; 
    
    # Examples
    ```jldoctest
    julia> ProtoSyn.Peptides.Calculators.Caterpillar.get_default_caterpillar_solvation_energy()
    🞧  Energy Function Component:
    +---------------------------------------------------+
    | Name           | Caterpillar_Solv                 |
    | Alpha (α)      | 1.0                              |
    | Update forces  | false                            |
    | Calculator     | neighbour_vector                 |
    +---------------------------------------------------+
     |    +----------------------------------------------------------------------------------+
     ├──  ● Settings                      | Value                                            |
     |    +----------------------------------------------------------------------------------+
     |    | hydrophobicity_weight         | nv_scalling_all_contributions                    |
     |    | identification_curve          | sigmoid                                          |
     |    | Ω                             | 24.0                                             |
     |    | sc                            | 1.0                                              |
     |    | hydrophobicity_map            | Dict{String, Float64}(21 components)             |
     |    | rmax                          | 9.0                                              |
     |    +----------------------------------------------------------------------------------+
     |    
     └──  ●  Selection:
          └── BinarySelection ❯  | "or" (Atom)
               ├── FieldSelection › Atom.name = CB
               └── BinarySelection ❯  & "and" (Atom)
                    ├── FieldSelection › Atom.name = HA2
                    └── FieldSelection › Residue.name = GLY
    ```
    """
    function get_default_caterpillar_solvation_energy(;α::T = 1.0) where {T <: AbstractFloat}
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