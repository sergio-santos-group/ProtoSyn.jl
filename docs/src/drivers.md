# Drivers

Drivers are a set set of functions that drive the simulation to new states, often exploring the conformational space of the system in Monte Carlo algorithms (with the combination of [Mutators](@ref Mutators) with acceptance/rejection rules) or the application of minimizers (such as the Steppest Descent Algorithm).

## Monte Carlo

```@meta
CurrentModule = Drivers.MonteCarlo
```

This section provides a description on the Monte Carlo Driver. This Driver iterates over a set amount of `n_steps` (defined in the [`MonteCarloDriver`](@ref)), sampling new conformations to the [`Common.State`](@ref) and accepting or rejecting them based on the Metropolis Algorithm.

```@docs
MonteCarloDriver
run!
```

## Steepest Descent

```@meta
CurrentModule = Drivers.SteepestDescent
```

This section provides a description on the Steepest Descent Driver. This Driver attempts to minimize the system energy based on the provided [`Topology`](@ref Forcefield), as it calculates the forces acting on each atom according to the defined [`SteepestDescentDriver`](@ref). The finalization criteria is:
- Maximum number of steps was achieved (SteepestDescentDriver.n_steps).
- Maximum force calculated is below the force tolerance (SteepestDescentDriver.f_tol).
- The gamma (γ) applied in the next step of the Steepest Descent is below machine precision.

```@docs
SteepestDescentDriver
run!
```