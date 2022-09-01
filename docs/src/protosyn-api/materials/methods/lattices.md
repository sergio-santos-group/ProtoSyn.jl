```@meta
CurrentModule = ProtoSyn.Materials
```

# Lattices

In the next section, a brief description of the available methods to generate new atomic lattices is provided. These lattices can be useful, among other things, to benchmark and test new [`EnergyFunction`](@ref ProtoSyn.Calculators.EnergyFunction) and [`EnergyFunctionComponent`](@ref ProtoSyn.Calculators.EnergyFunctionComponent) instances, for example. 

```@docs
primitive
```

![ProtoSyn Primitive Lattice](../../../assets/ProtoSyn-primitive-lattice.gif)

**Figure 1 |** An example of the [`primitive`](@ref) lattice.

```@docs
body_centered
face_centered
```

![ProtoSyn Body Centered Lattice](../../../assets/ProtoSyn-body-centered-lattice.gif)

**Figure 2 |** An example of the [`body_centered`](@ref) lattice.