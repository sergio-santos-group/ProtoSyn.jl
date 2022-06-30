```@meta
CurrentModule = ProtoSyn
```

# Common

The [Common](@ref) module is a special type of module, in ProtoSyn. Its objective is to group commonly performed tasks, such as defining a default [`EnergyFunction`](@ref ProtoSyn.Calculators.EnergyFunction) instances or defining commonly used [`Callback`](@ref ProtoSyn.Drivers.Callback) instances.

## Default Energy Function

```@docs
Common.default_energy_function
```

## Default Callbacks

```@docs
Common.default_energy_step_callback
Common.default_energy_step_frame_callback
Common.default_energy_step_detailed
Common.default_energy_step_frame_detailed
Common.detailed_message
```