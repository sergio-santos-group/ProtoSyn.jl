```@meta
CurrentModule = ProtoSyn.Calculators.Restraints
```

# Bond Distance Restraint

The [Bond Distance Restraint](@ref) is a potential restraint. However, in contrast with other [Potential Restraints](@ref calculators-potential-restraints) in ProtoSyn, a [`distance_matrix`](@ref ProtoSyn.Calculators.distance_matrix) is not calculated. Instead, the [`Pose`](@ref)'s [Graph](@ref graph-types) is iterated (See [Counters and Iterators](@ref)), and all bonds of an [`Atom`](@ref) are measured and evaluated according to a given potential.

```@docs
calc_bond_distance_restraint
get_default_bond_distance_restraint
```

![ProtoSyn Bond Distance Restraint](../../../assets/ProtoSyn-bond-distance-restraint.png)

**Figure 1 |** A diagram representation of the [Bond Distance Restraint](@ref) [`EnergyFunctionComponent`](@ref). The applied potential is a quadratic function, with domain between `[x0, +∞]`. The `x0` setting is parametrized based on observed distances between paits of [`Atom`](@ref) instances, according to the provided table (saved in `ProtoSyn.Units.max_bond_lengths`). If any given pair of atoms is not found in this table, the default `x0` value of the [`EnergyFunctionComponent`](@ref) is used instead (2.0 Å, by default).