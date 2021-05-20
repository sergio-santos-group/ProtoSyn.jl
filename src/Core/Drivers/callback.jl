using Printf

"""
    Callback(event::Function, frequency::Int)

A [`Callback`](@ref) instance is a _functor_, called with the following
signature:

```
(callback::Callback)(pose::Pose, driver_state::ProtoSyn.Drivers.DriverState)
```

A [`Callback`](@ref) expects the given `DriverState` `driver_state` to have a
`:step` field. The input [`Pose`](@ref) `pose` and `DriverState` `driver_state`
are passed as arguments to the `event` method every `callback.frequency` steps.

# Fields
* `event::Function` - Receives a [`Pose`](@ref) `pose` and `DriverState` `driver_state` as input arguments, should print some useful information to the user;
* `frequency::Int` - the `:event` method is called every `:frequency` steps (when comapred to the given `DriverState` `driver_state.step`).

# Examples
```jldoctest
julia> cb = ProtoSyn.Drivers.Callback((pose::Pose, driver_state::ProtoSyn.Drivers.DriverState) -> println("Status OK"), 1)
+----------------------------------------------------------+
| Index | Field           | Value                          |
+----------------------------------------------------------+
| 1     | Frequency       | 1                              |
| 2     | Event           | #5                             |
+----------------------------------------------------------+

julia> cb(pose, driver_state)
Status OK
```
"""
mutable struct Callback
    event::Function # Receives a pose and a DriverState
    frequency::Int
end

function (callback::Callback)(pose::Pose, driver_state::ProtoSyn.Drivers.DriverState)
    if driver_state.step % callback.frequency == 0
        callback.event(pose, driver_state)
    end
end

function Base.show(io::IO, cb::Callback)
    println(io, "+"*repeat("-", 58)*"+")
    @printf(io, "| %-5s | %-15s | %-30s |\n", "Index", "Field", "Value")
    println(io, "+"*repeat("-", 58)*"+")
    @printf(io, "| %-5d | %-15s | %-30d |\n", 1, "Frequency", cb.frequency)
    @printf(io, "| %-5d | %-15s | %-30s |\n", 2, "Event", "$(cb.event)")
    println(io, "+"*repeat("-", 58)*"+")
end