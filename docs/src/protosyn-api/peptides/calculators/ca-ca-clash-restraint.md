```@meta
CurrentModule = ProtoSyn.Peptides.Calculators.Restraints
```

# Cα - Cα Clash Restraint

The [Cα - Cα Clash Restraint](@ref) is a potential restraint (See [Potential restraints](@ref calculators-potential-restraints) for more details), applying a penalty for pairs of interacting Cα atoms with distances lower than a set value, thus reducing the occurrence of [`Residue`](@ref) overlap.

```@docs
get_default_ca_clash_restraint
```

![ProtoSyn CA-CA Restraint](../../../assets/ProtoSyn-ca-ca-clash-restraint.png)

**Figure 1 |** A diagram representation of the [Cα - Cα Clash Restraint](@ref) [`EnergyFunctionComponent`](@ref ProtoSyn.Calculators.EnergyFunctionComponent). The applied potential is a flat bottom function, with only the lower bound having and energy penatly. The `d1` and `d2` settings can be customly parametrized. This energy is only calculated for pair of interacting Cα [`Atom`](@ref) instances (in different [`Residue`](@ref) instances, as defined the the provided `AbstractSelection` and mask).