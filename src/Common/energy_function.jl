# COMMON ENERGY FUNCTIONS --------------------------------------------------

export default_energy_function

@doc """
    default_energy_function(::Type{T}) where {T <: AbstractFloat}
    default_energy_function()

Returns a default energy function for ProtoSyn. As of 2021, this includes
the following terms:
+----------------------------------------------------------+
| Index | Component name                      | Weight (ɑ) |
+----------------------------------------------------------+
| 1     | Caterpillar_Solvation               | 0.010      |
| 2     | TorchANI_ML_Model                   | 1.000      |
| 3     | Bond_Distance_Restraint             | 1.000      |
| 4     | Clash_Restraint                     | 100.000    |
+----------------------------------------------------------+
**Note:** If no Type{T} is provided, will use ProtoSyn.Units.defaultFloat;

# Examples
```jldoctest
julia> ProtoSyn.Common.default_energy_function()
...
```
"""
default_energy_function(::Type{T}) where {T <: AbstractFloat} = begin
    return Calculators.EnergyFunction(Dict(
        Calculators.TorchANI.torchani_model => T(1.0),
        Peptides.Calculators.Caterpillar.solvation_energy => T(0.01),
        Calculators.Restraints.bond_distance_restraint => T(1.0),
        ProtoSyn.Calculators.Restraints.clash_restraint => T(100.0)
    ))
end

default_energy_function() = begin
    default_energy_function(defaultFloat)
end

# ------------------------------------------------------------------------------

export xmlrpc_energy_function

@doc """
    xmlrpc_energy_function(::Type{T}) where {T <: AbstractFloat}
    xmlrpc_energy_function()

Returns a default energy function for ProtoSyn. As of 2021, this includes
the following terms:
+----------------------------------------------------------+
| Index | Component name                      | Weight (ɑ) |
+----------------------------------------------------------+
| 1     | TorchANI_ML_Model_XMLRPC            | 1.000      |
| 2     | Clash_Restraint                     | 100.000    |
| 3     | Bond_Distance_Restraint             | 1.000      |
| 4     | Caterpillar_Solvation               | 0.010      |
+----------------------------------------------------------+
**Note:** If no Type{T} is provided, will use ProtoSyn.Units.defaultFloat;
**Note:** This energy function employs the XML-RPC protocol whenever possible
(is slower, but safe in term of CUDA running out of memory).

# Examples
```jldoctest
julia> ProtoSyn.Common.xmlrpc_energy_function()
...
```

# See also
`default_energy_function`
"""
xmlrpc_energy_function(::Type{T}) where {T <: AbstractFloat} = begin
    return Calculators.EnergyFunction(Dict(
        Calculators.TorchANI.torchani_model_xmlrpc => T(1.0),
        Peptides.Calculators.Caterpillar.solvation_energy => T(0.01),
        Calculators.Restraints.bond_distance_restraint => T(1.0),
        ProtoSyn.Calculators.Restraints.clash_restraint => T(100.0)
    ))
end

xmlrpc_energy_function() = begin
    xmlrpc_energy_function(defaultFloat)
end