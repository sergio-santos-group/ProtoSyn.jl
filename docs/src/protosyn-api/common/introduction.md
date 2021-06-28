```@meta
CurrentModule = ProtoSyn.Common
```

# Common

The [Common](@ref) module is a pecial type of module, in ProtoSyn. Its function is to group commonly performed tasks, such as defining a default [`EnergyFunction`](@ref ProtoSyn.Calculators.EnergyFunction) instances or defining commonly used [`Callback`](@ref ProtoSyn.Drivers.Callback) instances.

## Default Energy Function

```@docs
default_energy_function
default_xmlrpc_energy_function
```

## Default Callbacks

```@docs
default_energy_step_callback
default_energy_step_frame_callback
```