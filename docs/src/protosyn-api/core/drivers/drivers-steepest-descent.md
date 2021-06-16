```@meta
CurrentModule = ProtoSyn.Drivers
```

# Steepest Descent Driver

The [Steepest Descent Driver](@ref) introduces a simple gradient descent, calculating the forces acting on each atom and applying a relaxing translational movement, until convergence or end of simulation.

```@docs
SteepestDescent
SteepestDescentState
```

![ProtoSyn steepest-descent-driver](../../../assets/ProtoSyn-steepest-descent-driver.png)

**Figure 1 |** A schematic overview of the different stages in a [`SteepestDescent`](@ref) `Driver`: **(1)** - Calculate the forces being felt on the system; **(2)** - Apply the forces, multiplied by the current stepsize; **(3)** - Adjust the stepsize, based on the new system's energy.