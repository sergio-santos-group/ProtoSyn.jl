```@meta
CurrentModule = ProtoSyn
```

# Hydrogen bonds

ProtoSyn aims to evaluate the stabilizing energy contribution of hydrogen bonds using a simple geometric potential. For each [`HydrogenBondPair`](@ref Calculators.HydrogenBonds.HydrogenBondPair) in a pre-calculated [`HydrogenBondNetwork`](@ref Calculators.HydrogenBonds.HydrogenBondNetwork), the hydrogen bond energetic contribution is given by product of a given distance-based potential (see [Available potentials](@ref)) by the $cos(α)cos(β)$, where $α$ and $β$ are the angles between the donor and receptor [`Atom`](@ref) instances (in each "side" of the interaction).

![ProtoSyn Hydrogen Bonds](../../../assets/ProtoSyn-hydrogen-bonds.png)

**Figure 1 |** A schematic representation of the geometric criteria taken into account when estimating the hydrogen bonds stabilizing interactions. Besides the inter-atomic distance, ProtoSyn also calculates the $α$ and $β$ angles, with maximum energetic reward being given to planar conformations of the dimeric system.

The [Hydrogen bonds](@ref) section is subdivided in the following sub-sections, for organizational purposes:

+ [Hydrogen bond network prediction](@ref)
+ [Hydrogen bonds EnergyFunctionComponent](@ref)

## Hydrogen bond network prediction

The [`generate_hydrogen_bond_network`](@ref Calculators.HydrogenBonds.generate_hydrogen_bond_network) attempts to identify all donor and receptor [`Atom`](@ref) instances in a given [`Pose`](@ref) using rudimentary criteria. This may change in future versions of ProtoSyn.

```@docs
Calculators.HydrogenBonds.HydrogenBondPair
Calculators.HydrogenBonds.HydrogenBondNetwork
Calculators.HydrogenBonds.generate_hydrogen_bond_network
```

## Hydrogen bonds EnergyFunctionComponent

By default, ProtoSyn generates a new [`HydrogenBondNetwork`](@ref Calculators.HydrogenBonds.HydrogenBondNetwork) each call to the default [Hydrogen bonds EnergyFunctionComponent](@ref). This may cause performance losses. If no mutation/design process is being performed, consider employing [`fixate_hydrogen_bond_network!`](@ref Calculators.HydrogenBonds.fixate_hydrogen_bond_network!) to set a static [`HydrogenBondNetwork`](@ref ProtoSyn.Calculators.HydrogenBonds.HydrogenBondNetwork).

```@docs
Calculators.HydrogenBonds.calc_hydrogen_bond_network
Calculators.HydrogenBonds.get_default_hydrogen_bond_network
Calculators.HydrogenBonds.fixate_hydrogen_bond_network!
```