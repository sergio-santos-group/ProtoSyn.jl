```@meta
CurrentModule = ProtoSyn.Peptides.Calculators.Restraints
```

# Sidechain Clash Restraint

The [Sidechain Clash Restraint](@ref) is a potential restraint (See [Potential restraints](@ref calculators-potential-restraints) for more details), applying a penalty for pairs of interacting atoms (belonging to different sidechains) with distances under a given threshold.

```@docs
get_default_sidechain_clash_restraint
```

![ProtoSyn Sidechain Clash Restraint](../../../assets/ProtoSyn-sidechain-clash-restraint.png)

**Figure 1 |** A diagram representation of the [Sidechain Clash Restraint](@ref) [`EnergyFunctionComponent`](@ref ProtoSyn.Calculators.EnergyFunctionComponent). The applied potential is a flat bottom function, with only the lower bound having and energy penatly. The `d1` and `d2` settings can be customly parametrized. This energy is only calculated for pair of interacting [`Atom`](@ref) instances in sidechains of different [`Residue`](@ref) instances (as defined the the provided `AbstractSelection` and mask).