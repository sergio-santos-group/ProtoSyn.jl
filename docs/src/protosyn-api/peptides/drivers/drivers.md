```@meta
CurrentModule = ProtoSyn.Peptides.Drivers
```

# Drivers

Each module in ProtoSyn may include a [Drivers section](@ref), adding new `Driver` types and methods. In the next sections, a closer look at the `Driver` instances introduced in the Peptides module is provided.

```@docs
RotamerBlitz
RotamerBlitzState
```

![ProtoSyn rotamer-blitz-driver](../../../assets/ProtoSyn-rotamer-blitz-driver.png)

**Figure 1 |** A schematic overview of the different stages in a [`RotamerBlitz`](@ref) `Driver`: **(1)** - For a single [`Residue`](@ref), all N most likely [`Rotamer`](@ref ProtoSyn.Peptides.Rotamer) instances are applied (using the [`apply!`](@ref ProtoSyn.Peptides.apply!) method) and the resulting conformation's energy is evaluated; **(2)** - After trying all [`Rotamer`](@ref ProtoSyn.Peptides.Rotamer) instances, the less energetic conformation is recovered and applied. This is performed to all [`Residue`](@ref) instances in a given [`Pose`](@ref), in random order; **(3)** - This process continues for a set amount of steps, resulting in a lower energy conformation for the sidechains in a [`Pose`](@ref).