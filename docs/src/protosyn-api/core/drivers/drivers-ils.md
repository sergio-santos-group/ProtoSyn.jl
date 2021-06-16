```@meta
CurrentModule = ProtoSyn.Drivers
```

# ILS Driver

The [ILS Driver](@ref) performs the Iterated Local Search algorithm, introducing a large conformational change in between local optimization runs.

```@docs
ILS
ILSState
```

![ProtoSyn ils-driver](../../../assets/ProtoSyn-ils-driver.png)

**Figure 1 |** A schematic overview of the different stages in a [`ILS`](@ref) `Driver`: **(1)** - A `jump!` introduces a relative large conformational change into the system's structure, taking the [`Pose`](@ref) into a new conformational space; **(2)** - An `inner_driver!` optimizes the new structure, finding new local minima.