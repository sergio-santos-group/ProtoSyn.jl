```@meta
CurrentModule = ProtoSyn.Peptides.Calculators.Caterpillar
```

# Caterpillar Solvation

The [Caterpillar Solvation](@ref) energy is a coarse-grained approximation to a solvation energy energetic contribution, introducing both a penalty for the exposure of hydrophobic [`Residue`](@ref) instances and a penalty for the burial of hydrophylic aminoacids. 

```@docs
calc_solvation_energy
get_default_solvation_energy
```

![ProtoSyn Caterpillar Solvation](../../../assets/ProtoSyn-caterpillar-solvation.png)

**Figure 1 |** A diagram representation of the [Caterpillar Solvation](@ref) [`EnergyFunctionComponent`](@ref ProtoSyn.Calculators.EnergyFunctionComponent). The energy penalty is proportional to the hydrophobicity value (in the hydrophobicity map) multiplied by the excess number of Cα contacts (above Ω).