```@meta
CurrentModule = ProtoSyn
```

# [SASA](@id calculators-sasa)

SASA stands for Solvent Accessible Surface Area and is a measure of the total surface area of a given solute. Multiple approaches have been proposed for the accurate estimation of this property. ProtoSyn, by default, employs the Overlapping Spheres (OLS) algorithm, a variant of the Shrake and Rupley algorithm (for more details, see [https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2712621/](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2712621/) and the [`calc_sasa`](@ref Calculators.SASA.calc_sasa) docstring). This section includes the following topics:

+ [SASA EnergyFunctionComponent](@ref)
+ [SASA solvation EnergyFunctionComponent](@ref)

## SASA EnergyFunctionComponent

The following [`EnergyFunctionComponent`](@ref ProtoSyn.Calculators.EnergyFunctionComponent) returns an energy value proportional to the system's SASA.

```@docs
Calculators.SASA.calc_sasa
Calculators.SASA.get_default_sasa
```

!!! ukw "Note:"
    Modern implicit solvation models often employ a hybrid approach, dubbed "SASA/GB". In short, the Generalized Born model attempts to estimate the enthalpic contribution of solvating a given molecule, and the SASA model calculates the entropy contribution of "opening space" for the solvation of such a molecule. Both models act together to provide a more clear picture of the solvation potential of that system. Consider employing the [`get_default_sasa`](@ref ProtoSyn.Calculators.SASA.get_default_sasa) [`EnergyFunctionComponent`](@ref ProtoSyn.Calculators.EnergyFunctionComponent) in conjunction with the [`get_default_gb`](@ref ProtoSyn.Calculators.GB.get_default_gb) [`EnergyFunctionComponent`](@ref ProtoSyn.Calculators.EnergyFunctionComponent).

## SASA solvation EnergyFunctionComponent

In contrast with the [`get_default_sasa`](@ref ProtoSyn.Calculators.SASA.get_default_sasa) [`EnergyFunctionComponent`](@ref ProtoSyn.Calculators.EnergyFunctionComponent), the [SASA solvation EnergyFunctionComponent](@ref) attempts to attribute a solvation energy value to each [`Residue`](@ref) in the system, based on an hydrophobicity map. In sum, exposed hydrophobic [`Residue`](@ref) instances are penalized while exposed hydrophilic [`Residue`](@ref) instances are rewarded (and vice-versa for buried [`Residue`](@ref) instances).

```@docs
Calculators.SASA.calc_sasa_energy
Calculators.SASA.get_default_sasa_energy
```