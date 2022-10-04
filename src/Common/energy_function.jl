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
| 3     | Bond_Distance_Rest                            |       1.00   |
| 4     | Hydrogen_Bonds                                |       0.50   |
| 5     | Coulomb                                       |       0.60   |
| 6     | GB_Solvation                                  |       1.00   |
| 7     | SASA                                          |   5.00e-04   |
+----------------------------------------------------------------------+
 ● Update forces: false
 ● Selection: Set
 └── TrueSelection (Atom)
```
"""
default_energy_function(::Type{T}) where {T <: AbstractFloat} = begin
    torchani    = ProtoSyn.Calculators.TorchANI.get_default_torchani_model(α = 0.0003)
    gb          = ProtoSyn.Calculators.GB.get_default_gb(α = 0.5)
    sasa        = ProtoSyn.Peptides.Calculators.SASA.get_default_sasa(α = 0.003)
    coulomb     = ProtoSyn.Calculators.Electrostatics.get_default_coulomb(α = 0.6)
    hydro_bonds = ProtoSyn.Calculators.HydrogenBonds.get_default_hydrogen_bond_network(α = 0.5)
    atom_clash  = ProtoSyn.Calculators.Restraints.get_default_all_atom_clash_restraint(α = 1.0)
    bond        = ProtoSyn.Calculators.Restraints.get_default_bond_distance_restraint(α = 0.05)
    
    return ProtoSyn.Calculators.EnergyFunction([
        torchani,
        atom_clash,
        bond,
        hydro_bonds,
        coulomb,
        gb,
        sasa])
end

default_energy_function() = begin
    default_energy_function(ProtoSyn.Units.defaultFloat)
end