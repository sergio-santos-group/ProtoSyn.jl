```@meta
CurrentModule = ProtoSyn.Peptides.Calculators.Restraints
```

# Contact Restraint

The [Contact Restraint](@ref) is a potential restraint (See [Potential restraints](@ref calculators-potential-restraints) for more details), applying a penalty for pairs of interacting Cα atoms with distances larger than a set value defined in a contact map.

```@docs
get_default_contact_restraint
```

![ProtoSyn Contact Restraint](../../../assets/ProtoSyn-contact-restraint.png)

**Figure 1 |** A diagram representation of the [Contact Restraint](@ref) [`EnergyFunctionComponent`](@ref ProtoSyn.Calculators.EnergyFunctionComponent). The applied potential is a flat bottom function, with only the upper bound having and energy penatly. The `d3` and `d4` settings can be customly parametrized. This energy is only calculated for pair of interacting Cα [`Atom`](@ref) instances (in different [`Residue`](@ref) instances, as defined the the provided `AbstractSelection` and mask). The penalty scale is defined in the loaded contact map file.