# COMMON ENERGY FUNCTIONS --------------------------------------------------

export default_energy_function

@doc """
    default_energy_function(::Type{T}) where {T <: AbstractFloat}
    default_energy_function()

Returns a default [`EnergyFunction`](@ref ProtoSyn.Calculators.EnergyFunction)
for ProtoSyn. If no `Type{T}` is provided, will use
`ProtoSyn.Units.defaultFloat`. This
[`EnergyFunction`](@ref ProtoSyn.Calculators.EnergyFunction) is optimized for
evaluating protein systems.

# Examples
```jldoctest
julia> ProtoSyn.Common.default_energy_function()
🗲  Energy Function (7 components):
+----------------------------------------------------------------------+
| Index | Component name                                | Weight (α)   |
+----------------------------------------------------------------------+
| 1     | TorchANI_ML_Model                             |       1.00   |
| 2     | All_Atom_Clash_Rest                           |       1.00   |
| 3     | Hydrogen_Bonds                                |       0.50   |
| 4     | Coulomb                                       |       0.60   |
| 5     | TorchANI_Ref_Energy                           |       1.00   |
| 6     | GB_Solvation                                  |       1.00   |
| 7     | SASA_Solvation                                |   5.00e-04   |
+----------------------------------------------------------------------+
 ● Update forces: false
 ● Selection: Set
 └── TrueSelection (Atom)
```
"""
default_energy_function(::Type{T}) where {T <: AbstractFloat} = begin
    torchani = ProtoSyn.Calculators.TorchANI.get_default_torchani_model()

    ref = ProtoSyn.Calculators.TorchANI.get_default_torchani_internal_energy()
    ref.settings[:use_ensemble] = false
    
    gb = ProtoSyn.Calculators.GB.get_default_gb(α = 1.0)
    gb.selection            = SidechainSelection()
    gb.settings[:ϵ_protein] = 25.0
    
    sasa = ProtoSyn.Peptides.Calculators.SASA.get_default_sasa_energy(α = 0.0005)
    max_sasas = copy(ProtoSyn.Peptides.Calculators.SASA.default_sidechain_max_sasa)
    max_sasas["CBZ"] = 1217.0
    sasa.settings[:max_sasas] = max_sasas
    
    coulomb = ProtoSyn.Calculators.Electrostatics.get_default_coulomb(α = 0.6)
    coulomb.settings[:mask]      = ProtoSyn.Calculators.get_intra_residue_mask
    coulomb.settings[:potential] = ProtoSyn.Calculators.get_bump_potential_charges(c = 0.0, r = 20.0)
    
    hydro_bonds = ProtoSyn.Calculators.HydrogenBonds.get_default_hydrogen_bond_network(α = 0.5)
    
    atom_clash = ProtoSyn.Calculators.Restraints.get_default_all_atom_clash_restraint(α = 1.0)
    atom_clash.settings[:mask] = ProtoSyn.Calculators.get_intra_residue_mask
    
    return ProtoSyn.Calculators.EnergyFunction([
        torchani,
        atom_clash,
        hydro_bonds,
        coulomb,
        ref,
        gb,
        sasa])
end

default_energy_function() = begin
    default_energy_function(ProtoSyn.Units.defaultFloat)
end