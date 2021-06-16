# COMMON ENERGY FUNCTIONS --------------------------------------------------

export default_energy_function

@doc """
    default_energy_function(::Type{T}) where {T <: AbstractFloat}
    default_energy_function()

Returns a default energy function for ProtoSyn. As of ProtoSyn 0.4, this
includes the following terms:
```
+----------------------------------------------------------+
| Index | Component name                      | Weight (α) |
+----------------------------------------------------------+
| 1     | TorchANI_ML_Model                   | 1.000      |
| 2     | Caterpillar_Solvation               | 0.010      |
| 3     | Bond_Distance_Restraint             | 1.000      |
| 4     | Cα-Cα_Clash_Restraint               | 100.000    |
+----------------------------------------------------------+
```
!!! note 
    If no Type{T} is provided, will use ProtoSyn.Units.defaultFloat;

# Examples
```jldoctest
julia> ProtoSyn.Common.default_energy_function()
🗲  Energy Function (4 components):
+----------------------------------------------------------------------+
| Index | Component name                                | Weight (α)   |
+----------------------------------------------------------------------+
| 1     | TorchANI_ML_Model                             |      1.000   |
| 2     | Caterpillar_Solvation                         |      0.010   |
| 3     | Bond_Distance_Restraint                       |      1.000   |
| 4     | Cα-Cα_Clash_Restraint                         |    100.000   |
+----------------------------------------------------------------------+
```
"""
default_energy_function(::Type{T}) where {T <: AbstractFloat} = begin
    return ProtoSyn.Calculators.EnergyFunction([
        ProtoSyn.Calculators.TorchANI.get_default_torchani_model(α = 1.0),
        ProtoSyn.Peptides.Calculators.Caterpillar.get_default_solvation_energy(α = 0.01),
        ProtoSyn.Calculators.Restraints.get_default_bond_distance_restraint(α = 1.0),
        ProtoSyn.Peptides.Calculators.Restraints.get_default_ca_clash_restraint(α = 100.0)
    ])
end

default_energy_function() = begin
    default_energy_function(ProtoSyn.Units.defaultFloat)
end

# ------------------------------------------------------------------------------

export xmlrpc_energy_function

@doc """
    xmlrpc_energy_function(::Type{T}) where {T <: AbstractFloat}
    xmlrpc_energy_function()

Returns a default energy function for ProtoSyn. As of ProtoSyn 0.4, this
includes the following terms:
```
+----------------------------------------------------------+
| Index | Component name                      | Weight (α) |
+----------------------------------------------------------+
| 1     | TorchANI_ML_Model_XMLRPC            | 1.000      |
| 2     | Caterpillar_Solvation               | 0.010      |
| 3     | Bond_Distance_Restraint             | 1.000      |
| 4     | Cα-Cα_Clash_Restraint               | 100.000    |
+----------------------------------------------------------+
```
!!! note 
    If no Type{T} is provided, will use ProtoSyn.Units.defaultFloat;
!!! note
    This energy function employs the XML-RPC protocol whenever possible (is slower, but safe in term of CUDA running out of memory).

# Examples
```jldoctest
julia> ProtoSyn.Common.xmlrpc_energy_function()
🗲  Energy Function (4 components):
+----------------------------------------------------------------------+
| Index | Component name                                | Weight (α)   |
+----------------------------------------------------------------------+
| 1     | TorchANI_ML_Model_XMLRPC                      |      1.000   |
| 2     | Caterpillar_Solvation                         |      0.010   |
| 3     | Bond_Distance_Restraint                       |      1.000   |
| 4     | Cα-Cα_Clash_Restraint                         |    100.000   |
+----------------------------------------------------------------------+

```

# See also
`default_energy_function`
"""
xmlrpc_energy_function(::Type{T}) where {T <: AbstractFloat} = begin
    return Calculators.EnergyFunction([
        ProtoSyn.Calculators.TorchANI.get_default_torchani_model_xmlrpc(α = 1.0),
        ProtoSyn.Peptides.Calculators.Caterpillar.get_default_solvation_energy(α = 0.01),
        ProtoSyn.Calculators.Restraints.get_default_bond_distance_restraint(α = 1.0),
        ProtoSyn.Peptides.Calculators.Restraints.get_default_ca_clash_restraint(α = 100.0)
    ])
end

xmlrpc_energy_function() = begin
    xmlrpc_energy_function(defaultFloat)
end