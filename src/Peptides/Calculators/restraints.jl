module Restraints

    using ProtoSyn
    using ProtoSyn.Calculators: EnergyFunctionComponent

    # * Default Energy Components ----------------------------------------------

    """
        get_default_sidechain_clash_restraint(;[α::T = 1.0], [mask::Opt{ProtoSyn.Mask} = nothing]) where {T <: AbstractFloat}

    Return the default sidechain clash restraint
    [`EnergyFunctionComponent`](@ref ProtoSyn.Calculators.EnergyFunctionComponent).
    `α` sets the component weight (on an
    [`EnergyFunction`](@ref ProtoSyn.Calculators.EnergyFunction) instance). If a
    `mask` is provided, the component will apply that mask every calculation
    (fixed mask) - recommended, except for design efforts. Otherwise, the
    default `get_intra_residue_mask` function will be used, which calculates a
    new intra-residue mask every calculation.

    # Sidechain clash energy settings
    - :d1, :d2, :d3, :d4 -> set each of the distances defining a flat-bottom potential (in Angstrom Å);
    - :selection -> defines the atom selection limiting the considered atoms for the calculation; 
    - :mask -> defines the mask applied to the energy and forces result;
    
    # See also
    [`ProtoSyn.Calculators.Restraints.calc_flat_bottom_restraint`](@ref)
    [`ProtoSyn.Calculators.get_intra_residue_mask`](@ref)

    # Examples
    ```
    julia> ProtoSyn.Peptides.Calculators.Restraints.get_default_sidechain_clash_restraint()
         Name : Clash_Sidechain_Restraint
    Weight(α) : 1.0
     Settings :
                :d4 => Inf
                :selection => UnarySelection ❯ ! "not" (Atom)
                              └── FieldSelection › Atom.name = r"^CA\$|^N\$|^C\$|^H\$|^O\$"
                :d2 => 1.5
              :mask => _intra_residue_mask
                :d1 => 1.0
                :d3 => Inf
    ```
    """
    function get_default_sidechain_clash_restraint(;α::T = ProtoSyn.Units.defaultFloat(1.0), mask::Opt{ProtoSyn.Mask} = nothing) where {T <: AbstractFloat}
        # * Note: The default :d1 and :d2 distances were parametrized based on
        # * the 2A3D PDB structure.
        
        _sele = !SidechainSelection()
        if mask === nothing
            mask = ProtoSyn.Calculators.get_intra_residue_mask
        end
        return EnergyFunctionComponent(
            "Sidechain_Clash",
            ProtoSyn.Calculators.Restraints.calc_flat_bottom_restraint,
            _sele,
            Dict{Symbol, Any}(:d1 => 1.0, :d2 => 1.5, :d3 => Inf, :d4 => Inf, :mask => mask),
            α,
            true)
    end


    """
        get_default_contact_restraint(filename::String; [α::T = 1.0]) where {T <: AbstractFloat}

    Return the default contact map restraint
    [`EnergyFunctionComponent`](@ref ProtoSyn.Calculators.EnergyFunctionComponent)
    by reading the given `filename`. `α` sets the component weight (on an
    [`EnergyFunction`](@ref ProtoSyn.Calculators.EnergyFunction)).

    !!! ukw "Note:"
        Since the contact map is fixed, any [`EnergyFunction`](@ref ProtoSyn.Calculators.EnergyFunction) containing this [`EnergyFunctionComponent`](@ref ProtoSyn.Calculators.EnergyFunctionComponent) can only be applied to one protein/sequence (no design). The attached contact map can be re-defined in `component.settings[:mask]`.
    
    !!! ukw "Note:"
        By default, this [`EnergyFunctionComponent`](@ref ProtoSyn.Calculators.EnergyFunctionComponent) does not calculate forces, as they would only be applied to the Cα [`Atom`](@ref) instances. This setting can be re-defined in `component.update_forces`.

    # Contact energy settings
    - :d1, :d2, :d3, :d4 -> set each of the distances defining a flat-bottom potential (in Angstrom Å);
    - :selection -> defines the atom selection limiting the considered atoms for the calculation; 
    - :mask -> defines the mask applied to the energy and forces result;
    
    # See also
    [`ProtoSyn.Calculators.Restraints.calc_flat_bottom_restraint`](@ref)
    [`ProtoSyn.Calculators.load_map`](@ref)

    # Examples
    ```
    julia> ProtoSyn.Peptides.Calculators.Restraints.get_default_contact_restraint("contact_map_example.txt")
    Name : Contact_Map
    Weight (α) : 1.0
 Update forces : false
       Setings :
            :d4 => 12.0
     :selection => FieldSelection › Atom.name = CA

            :d2 => 0.0
          :mask => Matrix{Float64}((73, 73)
            :d1 => 0.0
            :d3 => 8.0
    ```
    """
    function get_default_contact_restraint(;filename::String = "", α::T = ProtoSyn.Units.defaultFloat(1.0)) where {T <: AbstractFloat}
        _sele = an"CA"
        if filename === ""
            mask = nothing
        else
            mask  = ProtoSyn.Calculators.load_map(T, filename)
        end
        return EnergyFunctionComponent(
            "Contact_Map",
            ProtoSyn.Calculators.Restraints.calc_flat_bottom_restraint,
            _sele,
            Dict{Symbol, Any}(:d1 => 0.0, :d2 => 0.0, :d3 => 8.0, :d4 => 12.0, :mask => mask),
            α,
            false)
    end


    """
        get_default_ca_clash_restraint(;[α::T = 1.0]) where {T <: AbstractFloat}
    
    Return the default Cα-Cα clash restraint
    [`EnergyFunctionComponent`](@ref ProtoSyn.Calculators.EnergyFunctionComponent).
    `α` sets the component weight (on an
    [`EnergyFunction`](@ref ProtoSyn.Calculators.EnergyFunction)). If a `mask`
    is provided, the component will apply that mask every calculation (fixed
    mask) - recommended, except for design efforts. Otherwise, the default
    [`get_diagonal_mask`](@ref ProtoSyn.Calculators.get_diagonal_mask) function
    will be used, which calculates a new diagonal mask every calculation
    (effectly ignoring the same atom energetic contributions).

    !!! ukw "Note:"
        By default, this [`EnergyFunctionComponent`](@ref ProtoSyn.Calculators.EnergyFunctionComponent) does not calculate forces, as they would only be applied to the Cα atoms. This setting can be re-defined in `component.update_forces`.

    # Cα clash energy settings
    - :d1, :d2, :d3, :d4 -> set each of the distances defining a flat-bottom potential (in Angstrom Å);
    - :selection -> defines the atom selection limiting the considered atoms for the calculation; 
    - :mask -> defines the mask applied to the energy and forces result;
    
    # See also
    [`ProtoSyn.Calculators.Restraints.calc_flat_bottom_restraint`](@ref)
    [`ProtoSyn.Calculators.get_diagonal_mask`](@ref)

    # Examples
    ```
    julia> ProtoSyn.Peptides.Calculators.Restraints.get_default_ca_clash_restraint()
             Name : Cα-Cα_Clash_Restraint
       Weight (α) : 1.0
    Update forces : false
          Setings :
           :d4 => Inf
    :selection => FieldSelection › Atom.name = CA
           :d2 => 3.0
         :mask => _diagonal_mask
           :d1 => 1.0
           :d3 => Inf
    ```
    """
    function get_default_ca_clash_restraint(;α::T = ProtoSyn.Units.defaultFloat(1.0), mask::Opt{ProtoSyn.Mask} = nothing) where {T <: AbstractFloat}
        _sele = an"CA"
        if mask === nothing
            mask = ProtoSyn.Calculators.get_diagonal_mask
        end
        return EnergyFunctionComponent(
            "Cα-Cα_Clash_Restraint",
            ProtoSyn.Calculators.Restraints.calc_flat_bottom_restraint,
            _sele,
            Dict{Symbol, Any}(:d1 => 1.0, :d2 => 3.5, :d3 => Inf, :d4 => Inf, :mask => mask, :vlist => nothing),
            α,
            false)
    end
end