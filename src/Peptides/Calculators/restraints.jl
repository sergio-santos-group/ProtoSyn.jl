module Restraints

    using ProtoSyn
    using ProtoSyn.Calculators: EnergyFunctionComponent

    # * Default Energy Components ----------------------------------------------

    """
        get_default_sidechain_clash_restraint(;α::T = 1.0, mask::Opt{ProtoSyn.Mask} = nothing) where {T <: AbstractFloat}

    Return the default sidechain clash restraint `EnergyFunctionComponent`. `α`
    sets the component weight (on an `EnergyFunction`). If a `mask` is provided,
    the component will apply that mask every calculation (fixed mask) -
    recommended, except for design efforts. Otherwise, the default
    `get_intra_residue_mask` function will be used, which calculates a new
    intra-residue mask every calculation.

    # Sidechain clash energy settings
    - :d1, :d2, :d3, :d4 -> set each of the distances defining a flat-bottom potential (in Angstrom Å);
    - :selection -> defines the atom selection limiting the considered atoms for the calculation; 
    - :mask -> defines the mask applied to the energy and forces result;
    
    # See also
    `ProtoSyn.Calculators.Restraints.calc_flat_bottom_restraint`

    # Examples
    ```jldoctest
    julia> ProtoSyn.Peptides.Calculators.Restraints.get_default_sidechain_clash_restraint()
         Name : Clash_Sidechain_Restraint
    Weight(α) : 1.0
     Settings :
                :d4 => Inf
         :selection => UnarySelection{ProtoSyn.Stateless}(!, FieldSelection{ProtoSyn.Stateless,Atom}(r"^CA\$|^N\$|^C\$|^H\$|^O\$", :name, occursin))
                :d2 => 4.2
              :mask => _intra_residue_mask
                :d1 => 2.0
                :d3 => Inf
    ```
    """
    function get_default_sidechain_clash_restraint(;α::T = ProtoSyn.Units.defaultFloat(1.0), mask::Opt{ProtoSyn.Mask} = nothing) where {T <: AbstractFloat}
        # * Note: The default :d1 and :d2 distances were parametrized based on
        # * the 2A3D PDB structure.
        
        _sele = !an"^CA$|^N$|^C$|^H$|^O$"r
        if mask === nothing
            mask = ProtoSyn.Calculators.get_intra_residue_mask(_sele)
        end
        return EnergyFunctionComponent(
            "Clash_Sidechain_Restraint",
            ProtoSyn.Calculators.Restraints.calc_flat_bottom_restraint,
            Dict{Symbol, Any}(:d1 => 1.0, :d2 => 1.5, :d3 => Inf, :d4 => Inf, :selection => _sele, :mask => mask),
            α,
            true)
    end


    """
        get_default_contact_restraint(filename::String; α::T = 1.0) where {T <: AbstractFloat}

    Return the default contact map restraint `EnergyFunctionComponent` by
    reading the given `filename`. `α` sets the component weight (on an
    `EnergyFunction`). 
    !!! note
        Since the map is fixed, any energy function containing this component can only be applied to one protein/sequence. The attached contact map can be re-defined in `component.settings[:mask]`.
    !!! note
        By default, this component does not calculate forces, as they would only be applied to the Cα atoms. This setting can be re-defined in `component.update_forces`.

    # Contact energy settings
    - :d1, :d2, :d3, :d4 -> set each of the distances defining a flat-bottom potential (in Angstrom Å);
    - :selection -> defines the atom selection limiting the considered atoms for the calculation; 
    - :mask -> defines the mask applied to the energy and forces result;
    
    # See also
    `ProtoSyn.Calculators.Restraints.calc_flat_bottom_restraint`

    # Examples
    ```jldoctest
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
    function get_default_contact_restraint(filename::String; α::T = ProtoSyn.Units.defaultFloat(1.0)) where {T <: AbstractFloat}
        _sele = an"CA"
        mask  = ProtoSyn.Calculators.load_map(T, filename)
        return EnergyFunctionComponent(
            "Contact_Map",
            ProtoSyn.Calculators.Restraints.calc_flat_bottom_restraint,
            Dict{Symbol, Any}(:d1 => 0.0, :d2 => 0.0, :d3 => 8.0, :d4 => 12.0, :selection => _sele, :mask => mask),
            α,
            false)
    end

    """
        get_default_ca_clash_restraint(;α::T = ProtoSyn.Units.defaultFloat(1.0)) where {T <: AbstractFloat}
    
    Return the default Cα-Cα clash restraint `EnergyFunctionComponent`. `α`
    sets the component weight (on an `EnergyFunction`). If a `mask` is provided,
    the component will apply that mask every calculation (fixed mask) -
    recommended, except for design efforts. Otherwise, the default
    `get_diagonal_mask` function will be used, which calculates a new
    diagonal mask every calculation (effectly ignoring the same atom energetic
    contributions). 
    !!! note
        By default, this component does not calculate forces, as they would only be applied to the Cα atoms. This setting can be re-defined in `component.update_forces`.

    # Cα clash energy settings
    - :d1, :d2, :d3, :d4 -> set each of the distances defining a flat-bottom potential (in Angstrom Å);
    - :selection -> defines the atom selection limiting the considered atoms for the calculation; 
    - :mask -> defines the mask applied to the energy and forces result;
    
    # See also
    `ProtoSyn.Calculators.Restraints.calc_flat_bottom_restraint`

    # Examples
    ```jldoctest
    julia> ProtoSyn.Peptides.Calculators.Restraints.get_default_ca_clash_restraint()
             Name : Cα-Cα_Clash_Restraint
       Weight (α) : 1.0
    Update forces : false
          Setings :
           :d4 => Inf
    :selection => FieldSelection{ProtoSyn.Stateless,Atom}("CA", :name, isequal)
           :d2 => 3.0
         :mask => _diagonal_mask
           :d1 => 1.0
           :d3 => Inf
    ```
    """
    function get_default_ca_clash_restraint(;α::T = ProtoSyn.Units.defaultFloat(1.0), mask::Opt{ProtoSyn.Mask} = nothing) where {T <: AbstractFloat}
        _sele = an"CA"
        if mask === nothing
            mask = ProtoSyn.Calculators.get_diagonal_mask(_sele)
        end
        return EnergyFunctionComponent(
            "Cα-Cα_Clash_Restraint",
            ProtoSyn.Calculators.Restraints.calc_flat_bottom_restraint,
            Dict{Symbol, Any}(:d1 => 1.0, :d2 => 3.0, :d3 => Inf, :d4 => Inf, :selection => _sele, :mask => mask),
            α,
            false)
    end
end