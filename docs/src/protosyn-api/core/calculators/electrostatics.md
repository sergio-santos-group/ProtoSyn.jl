```@meta
CurrentModule = ProtoSyn
```

# [Electrostatics](@id calculators-electrostatics)

The following section goes into detail on how ProtoSyn aims to evaluate electrostatic interactions, divided in the following topics:

+ [Adding charges](@ref)
+ [Coulomb EnergyFunctionComponent](@ref)

## Adding charges

The first step in any eletrostatic energy calculation is the attribution of charges to each point-like particle (i.e.: each [`Atom`](@ref) instance). Users are free to manually define each [`Atom`](@ref)'s charge (by setting the corresponting `pose.state[atom].δ` field). More automation purposes, ProtoSyn makes available several methods to automatically set-up charges in a [`Pose`](@ref):

```@docs
Calculators.Electrostatics.assign_default_charges!
Calculators.Electrostatics.assign_acc2_eem_charges_from_file!
```

## Coulomb EnergyFunctionComponent

As long as a [`Pose`](@ref) has charges attributed, the default [Electrostatics](@ref calculators-electrostatics) coulomb [`EnergyFunctionComponent`](@ref ProtoSyn.Calculators.EnergyFunctionComponent) can be applied.

```@docs
Calculators.Electrostatics.calc_coulomb
Calculators.Electrostatics.get_default_coulomb
```