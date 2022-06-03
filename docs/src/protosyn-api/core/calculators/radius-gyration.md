```@meta
CurrentModule = ProtoSyn.Calculators
```

# Radius of gyration

The radius of gyration of a molecule is, in short, a measure of how "elongated" that structure is in each of the 3 dimensions (x, y and z). This can be useful, for example, in globular systems where an elongation in any of the dimensions is usually not likely to be a stabilizing conformational change. The following methods allow ProtoSyn to calculate the radius of gyration of a [`Pose`](@ref) and attribute a tentative energy value.

```@docs
Calculators.RG.calc_radius_gyration
Calculators.RG.calc_radius_gyration_energy
Calculators.RG.get_default_rg
```