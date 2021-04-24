```@meta
CurrentModule = ProtoSyn.Calculators
```

# Energy Function

This module contains types and methods that allow ProtoSyn to measure the fitness of a given [`State`](@ref), according to a specific [`EnergyFunction`](@ref) instance. The resulting energy value is stored in the [`State`](@ref)`.e` field, as a Dictionary of [`EnergyFunctionComponent`](@ref) instances. The `:total` energy of the system is the sum of each component contribution. In the next section, a further exploration on how to develop and set-up custom [`EnergyFunction`](@ref) instances is provided, organized in the following sections:

+ [Energy functions & energy function components](@ref)
+ [Creating a custom EnergyFunctionComponent](@ref)
+ [ProtoSyn acceleration types](@ref)
+ [Verlet lists](@ref)


### Energy functions & energy function components

An [`EnergyFunction`](@ref) is the main high-level object of ProtoSyn for energy calculation of a [`Pose`](@ref) instance, whose building blocks are the [`EnergyFunctionComponent`](@ref) instances.

```@docs
Calculators.EnergyFunctionComponent
Calculators.EnergyFunction
```

An [`EnergyFunction`](@ref) instance is actually a __functor__, meaning it is callable in order to apply any energy function type to any [`Pose`](@ref) instance, using the following signature:

```
(energy_function::EnergyFunction)(pose::Pose, update_forces::Bool = false)
```

An application example is provided bellow:

```@setup calculators
using ProtoSyn
using ProtoSyn.Peptides
using ProtoSyn.Calculators
res_lib = Peptides.grammar(Float64)
pose = ProtoSyn.build(res_lib, seq"GME")
torchani = Calculators.TorchANI.get_default_torchani_model()
energy_function = Calculators.EnergyFunction([torchani])
energy_function(pose)
```

```@repl calculators
torchani = Calculators.TorchANI.get_default_torchani_model()
energy_function = Calculators.EnergyFunction([torchani])
energy_function(pose)
energy_function(pose, true) # Calculating forces ...
pose.state.f
```

![ProtoSyn Energy Function](../../../assets/ProtoSyn-energy-function.png)

**Figure 1 |** A diagram representation of the [`EnergyFunction`](@ref) organization as a set of [`EnergyFunctionComponent`](@ref) instances. Each component contribution is summed to evaluate the `:total` energy of the system, according to each component's scaling factor `ɑ` (in this example, `1.0` for the *TorchANI Ensemble* [`EnergyFunctionComponent`](@ref) and `0.8` for the *Contact Map* [`EnergyFunctionComponent`](@ref)). Some components may be able to calculate the set of forces acting on each atom of the system, in which case this contribution can be toggled *on* and *off* by the `:update_forces` flag. Furthermore, each component may be parametrized by a set of options. In this example, the *Contact Map* [`EnergyFunctionComponent`](@ref) can be further parametrized by setting the `:d1`, `:d2`, `:d3` and `:d4` values (controlling the flat bottom restraint potential), as well as the `:selection ` on which to act and the `:mask` applied. Such settings are specific to each type of [`EnergyFunctionComponent`](@ref) and are explored in detail in each entry.

### Creating a custom EnergyFunctionComponent

As stated before, an [`EnergyFunctionComponent`](@ref) is mostly a smart wrap over a `calc` `Function` responsible for calculating a given energetic contribution based on a specific method, model or restraint type (as an example). Therefore, when developing custom [`EnergyFunctionComponent`](@ref) instances, certain aspects need to be taken into consideration. The following section lists particular restraints when creating custom [`EnergyFunctionComponent`](@ref) instances.

1. The `calc` `Function` signature:

When being called from an [`EnergyFunction`](@ref) instance, ProtoSyn expects the following `calc` `Function` signature.

```
calc(pose::Pose, update_forces::Bool = false; kwarg1::Type{kwarg1} = default_value1, kwarg2::Type{kwarg2} = default_value2, ...)
```

The `pose::Pose` and `update_forces::Bool = false` arguments are necessary and expected. Optionally, these can be followed by a variable number of arguments (kwargs), parameterizing and fine-tuning the usage of the `calc` `Function`. When calling this [`EnergyFunctionComponent`](@ref) from an [`EnergyFunction`](@ref),
these kwargs are passed from the [`EnergyFunctionComponent`](@ref)`.settings` dictionary, and can therefore be automatically and continuously adjusted in a simulation. Since it may happen that one or more kwargs may not be defined in the [`EnergyFunctionComponent`](@ref)`.settings`, defining a default value for each of the kwargs is considered a good practice when creating a custom [`EnergyFunctionComponent`](@ref).

Examples:

```
calc_torchani_ensemble(pose::Pose, update_forces::Bool = false)

calc_flat_bottom_restraint(pose::Pose, update_forces::Bool; d1::T = 0.0, d2::T = 0.0, d3::T = Inf, d4::T = Inf, selection::Opt{AbstractSelection} = nothing, mask::MaskMap = nothing) where {T <: AbstractFloat}
```

2. The [`EnergyFunctionComponent`](@ref)`.settings` field

As previously stated, the [`EnergyFunctionComponent`](@ref)`.settings` are optional. If present, these are passed to the [`EnergyFunctionComponent`](@ref)`.calc` `Function` as kwargs, and therefore should match and have the same name (as a `Symbol` instead). This allows for easy access and modification of internal parameters of the calculation on the fly.

Examples:

```
Dict{Symbol, Any}(:d1 => 1.0, :d2 => 3.0, :d3 => Inf, :d4 => Inf, :selection => an"CA", :mask => Calculators.get_diagonal_mask(an"CA"))
```

### ProtoSyn acceleration types

Certain [`EnergyFunctionComponent`](@ref) instances in ProtoSyn require specific software and hardware on the machine performing the calculation, or at least have specific implementations of the code based on the availability or not os said features (for example, certain functions have specific implementations for CUDA and GPU acceleration). For this reason, ProtoSyn defines 3 levels of acceleration (3 different acceleration types, loosely based on [Flynn's Taxonomy](https://www.geeksforgeeks.org/computer-architecture-flynns-taxonomy/)):

+ `SISD_0` - [Single Instruction Single Data](https://en.wikipedia.org/wiki/SISD): uses CPU, should work on most modern processors;
+ `SIMD_1` - [Single Instruction Multiple Data](https://en.wikipedia.org/wiki/SIMD): uses CPU, parallelizes a single instruction in a vector of data, available in modern processors with [MMX](https://wikipedia.org/wiki/MMX), [SSE](https://wikipedia.org/wiki/Streaming_SIMD_Extensions) and [AVX](https://wikipedia.org/wiki/Advanced_Vector_Extensions)/[AVX-512](https://en.wikipedia.org/wiki/AVX-512) technologies, among others;
+ `CUDA_2` - [Compute Unified Device Architecture](https://en.wikipedia.org/wiki/CUDA): parallel computing platform and API enabling GPU acceleration of code. If CUDA is available, this is the default acceleration type set up in ProtoSyn.

As start-up, ProtoSyn tries to identify the available software and hardware available and sets the fastest acceleration type possible as `ProtoSyn.acceleration.active`. This can, however, be overridden:

```@repl calculators
ProtoSyn.acceleration.active
ProtoSyn.acceleration.active = ProtoSyn.SIMD_1
```

### Verlet lists