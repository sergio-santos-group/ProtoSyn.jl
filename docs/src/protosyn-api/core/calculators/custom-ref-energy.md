```@meta
CurrentModule = ProtoSyn.Calculators
```

# Custom reference energy

Similarly to the logic behind [TorchANI reference energy EnergyFunctionComponent](@ref), it is often useful to attribute a custom reference energy to a [`Residue`](@ref) type (by name). This can be, for example, to promote a positive bias towards a certain sub-set of aminoacids during mutation/design processes. The following methods allow the introduction of custom reference energies in ProtoSyn energy evaluations.

```@docs
Calculators.calc_custom_ref_energy
Calculators.get_default_custom_ref_energy
```