# Drivers

Drivers are a set set of functions that drive the simulation to new states, often exploring the conformational space of the system in Monte Carlo algorithms (with the combination of [Mutators](@ref) with acceptance/rejection rules) or the application of minimizers (such as the Steppest Descent Algorithm).
For each Driver a set of required `callback functions` is further explained, detailing the necessary signature.

## Monte Carlo

```@meta
CurrentModule = Drivers.MonteCarlo
```

This section provides a description on the Monte Carlo Driver. This Driver iterates over a set amount of `steps` (defined in the [`ConfigParameters`](@ref)), sampling new conformations to the [`Common.State`](@ref) and accepting or rejecting them based on the Metropolis Algorithm.

```@docs
ConfigParameters
load_parameters
run!
```

## Steepest Descent

```@meta
CurrentModule = Drivers.SteepestDescent
```

This section provides a description on the Steepest Descent Driver. This Driver attempts to minimize the system energy based on the provided [`Topology`](@ref Forcefield), as it calculates the forces acting on each atom according to the defined [`ConfigParameters`](@ref). The finalization criteria is:
- Maximum number of steps was achieved (ConfigParameters.n_steps).
- Maximum force calculated is below the force tolerance (ConfigParameters.f_tol).

```@docs
ConfigParameters
load_parameters
run!
```