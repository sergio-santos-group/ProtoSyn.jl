module SASA

    using ProtoSyn
    using ProtoSyn.Peptides
    using YAML
    using Statistics

    """
        calculate_max_sasa(selection::AbstractSelection, residue_list::Vector{Char}; grammar::LGrammar = Peptides.grammar, probe_radius::T = 1.4, n_points::Int = 100) where {T <: AbstractFloat}
    
    For each [`Residue`](@ref) type in `residue_list`, build a template from the
    provided [`LGrammar`](@ref) `grammar` (`Peptides.grammar`, by default) and
    evaluate the maximum SASA (the single [`Residue`](@ref) instance completly
    solvated by the implicit solvent). The the given `AbstractSelection`
    `selection`, `probe_radius` and `n_points` settings are used to calculate
    the SASA value (see [`calc_sasa`](@ref ProtoSyn.Calculators.SASA.calc_sasa)).
    Note that, when using the maximum SASA for comparison with the evaluated
    SASA in a complex and partially occluded system (such as a protein
    environment), the same settings should be employed (`AbstractSelection`,
    `probe_radius` and `n_points`). Returns a dictionary with a maximum SASA
    value for each queried [`Residue`](@ref) type.

    !!! ukw "Note:"
        This method is used to define `default_sidechain_max_sasa` using the
        default settings from [`get_default_sasa_energy`](@ref):
        * `AbstractSelection`: (SidechainSelection() & !as"H") | an"CA"
        * `probe_radius`: 2.8 Å
        * `n_points`: 100

    # See also
    [`calculate_average_sasa`](@ref)

    # Examples
    ```jldoctest
    julia> ProtoSyn.Peptides.Calculators.SASA.calculate_max_sasa(((SidechainSelection() & !as"H") | an"CA"), ['R', 'G'], probe_radius = 2.8)
    Dict{String, Float64} with 2 entries:
      "GLY" => 100.0
      "ARG" => 260.0
    ```  
    """
    function calculate_max_sasa(selection::AbstractSelection, residue_list::Vector{Char}; grammar::LGrammar = Peptides.grammar, probe_radius::T = 1.4, n_points::Int = 100) where {T <: AbstractFloat}

        max_sasa = Dict{String, T}()
        for residue in residue_list
            res = string(residue)
            pose = Peptides.build(grammar, [res])
            sasa = ProtoSyn.Calculators.SASA.calc_sasa(pose, selection, false,
                probe_radius = probe_radius,
                n_points = n_points)[1]
            for entry in ProtoSyn.alt_residue_names[ProtoSyn.one_2_three[residue]]
                max_sasa[entry] = sasa
            end
        end

        return max_sasa
    end


    """
        calculate_average_sasa(pose::Pose, selection::AbstractSelection, max_sasas::Dict{String, T}; probe_radius::T = 1.4, n_points::Int = 100) where {T <: AbstractFloat}
    
    Calculates the average SASA exposure of the given [`Pose`](@ref) `pose` (as
    a % in relation to the maximum SASA, as defined in `max_sasas` - see
    [`calculate_max_sasa`](@ref)). If an `AbstractSelection` `selection` is
    provided, only the selected [`Atom`](@ref) instances will be considered for
    SASA calculation (using the given `probe_radius` and `n_points`). Note that
    for comparison, make sure the `max_sasas` were calcualted using the same 
    settings (`selection`, `probe_radius` and `n_points`). If using an
    `AbstractSelection` that only selects a given subset of [`Residue`](@ref)
    instances (instead of a subset of [`Atom`](@ref) instances in all
    [`Residue`](@ref) instances), note that the SASA can only be occluded by all
    other selected [`Residue`](@ref) instances. This generally tens to return a
    much higher average SASA value.

    !!! ukw "Note:"
        In the default Peptides SASA [`EnergyFunctionComponent`](@ref ProtoSyn.Calculators.EnergyFunctionComponent) (see [`get_default_sasa_energy`](@ref)), a `Ω` value sets the average SASA value for a [`Residue`](@ref) instance in the evaluated [`Pose`](@ref). [`calculate_average_sasa`](@ref) can help setting this value.

    # Examples
    ```
    julia> ProtoSyn.Peptides.Calculators.SASA.calculate_average_sasa(pose, ((SidechainSelection() & !as"H") | an"CA"), Peptides.Calculators.SASA.default_sidechain_max_sasa, probe_radius = 2.8)
    0.49055441526123295
    ``` 
    """
    function calculate_average_sasa(pose::Pose, selection::AbstractSelection, max_sasas::Dict{String, T}; probe_radius::T = 1.4, n_points::Int = 100) where {T <: AbstractFloat}

        s, _, sasas = ProtoSyn.Calculators.SASA.calc_sasa(pose, selection, false,
            probe_radius = probe_radius,
            n_points = n_points)

        a_resi_ids = [a.container.id for a in selection(pose, gather = true)]
        
        exposures = Vector{T}()
        # for residue in eachresidue(pose.graph)
        for residue in ProtoSyn.promote(selection, Residue)(pose, gather = true)
            resi_atoms_indexes = findall(==(residue.id), a_resi_ids)
            resi_sasa          = sum(sasas[resi_atoms_indexes])
            exposure           = resi_sasa / max_sasas[residue.name]
            push!(exposures, exposure)
        end

        return mean(exposures)
    end

    const default_sidechain_max_sasa = Dict{String, ProtoSyn.Units.defaultFloat}(
        YAML.load_file(
            joinpath(ProtoSyn.Peptides.resource_dir, "Calculators", "aa-max-sasa.yml")))

    """
        get_default_sasa_energy(;[α::T = 1.0]) where {T <: AbstractFloat}

    Return the default SASA energy [`EnergyFunctionComponent`](@ref ProtoSyn.Calculators.EnergyFunctionComponent). `α` sets
    the component weight (on an
    [`EnergyFunction`](@ref ProtoSyn.Calculators.EnergyFunction) instance, `1.0`
    by default). This function employs
    [`calc_sasa`](@ref ProtoSyn.Calculators.SASA.calc_sasa) as the `:calc`
    function.

    !!! ukw "Note:"
        This [`EnergyFunctionComponent`](@ref ProtoSyn.Calculators.EnergyFunctionComponent) is similar to the [`ProtoSyn.Calculators.SASA.get_default_sasa_energy`](@ref) version, but adds the correct `selection` and `max_sasas` dictionary (specifically for peptides).

    # Settings
    * `n_points::Int` - The number of points to generate in each [`Atom`](@ref) sphere (higher number of points leads to higher accuracy, at the expense of performance);
    * `probe_radius::T` - The distance of each point in a generated sphere to the central [`Atom`](@ref) instance. Any point within `probe_radius` of any other atom is considered buried [`Residue`](@ref) name (where T <: AbstractFloat);
    * `hydrophobicity_map::Dict{String, T}` - A dictionary of hydrophobicity values for each [`Residue`](@ref) name, positive values indicate hydrophobicity and vice-versa (where T <: AbstractFloat);
    * `max_sasas::Dict{String, T}` - A dictionary of max_sasa values (SASA values for fully-solvated [`Residue`](@ref) instances) for each [`Residue`](@ref) name (where T <: AbstractFloat);
    * `Ω::T` - The average exposure value (between 0.0 and 1.0), any SASA value bellow this percentage of max_sasa is considered buried (where T <: AbstractFloat) (see [`calculate_average_sasa`](@ref)).

    # See also
    [`calc_sasa_energy`](@ref ProtoSyn.Calculators.SASA.calc_sasa_energy)

    # Examples
    ```jldoctest
    julia> ProtoSyn.Peptides.Calculators.SASA.get_default_sasa_energy()
    🞧  Energy Function Component:
    +---------------------------------------------------+
    | Name           | SASA_Solvation                   |
    | Alpha (α)      | 1.0                              |
    | Update forces  | false                            |
    | Calculator     | calc_sasa_energy                 |
    +---------------------------------------------------+
     |    +----------------------------------------------------------------------------------+
     ├──  ● Settings                      | Value                                            |
     |    +----------------------------------------------------------------------------------+
     |    | max_sasas                     | Dict{String, Float64}(22 components)             |
     |    | hydrophobicity_map            | Dict{String, Float64}(22 components)             |
     |    | Ω                             | 0.5                                              |
     |    | probe_radius                  | 2.8                                              |
     |    | n_points                      | 100                                              |
     |    +----------------------------------------------------------------------------------+
     |    
     └──  ●  Selection:
          └── BinarySelection ❯  | "or" (Atom)
               ├── BinarySelection ❯  & "and" (Atom)
               |    ├── SidechainSelection › (Atom)
               |    └── UnarySelection ❯ ! "not" (Atom)
               |         └── FieldSelection › Atom.symbol = H
               └── FieldSelection › Atom.name = CA
    ```
    """
    function get_default_sasa_energy(;α::T = 1.0) where {T <: AbstractFloat}
        return ProtoSyn.Calculators.EnergyFunctionComponent(
            "SASA_Solvation",
            ProtoSyn.Calculators.SASA.calc_sasa_energy,
            (SidechainSelection() & !as"H") | an"CA",
            Dict{Symbol, Any}(
                :probe_radius                 => 2.8,
                :n_points                     => 100,
                :hydrophobicity_map           => Peptides.doolitle_hydrophobicity,
                :max_sasas                    => default_sidechain_max_sasa,
                :Ω                            => 0.5
            ),
            α,
            false
        )
    end


    """
    # TODO
    """
    function get_default_sasa(;α::T = 1.0) where {T <: AbstractFloat}
        return ProtoSyn.Calculators.EnergyFunctionComponent(
            "SASA",
            ProtoSyn.Calculators.SASA.calc_sasa,
            an"CA",
            Dict{Symbol, Any}(
                :n_points           => 100,
                :probe_radius       => 5.0,
            ),
            α,
            false)
    end
end