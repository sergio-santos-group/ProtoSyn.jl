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
✉  Callback:
+----------------------------------------------------------------------+
| Index | Field                     | Value                            |
+----------------------------------------------------------------------+
| 1     | Event                     | #1                               |
| 2     | Frequency                 | 1                                |
+----------------------------------------------------------------------+

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

function Base.show(io::IO, cb::Callback, level_code::Opt{LevelCode} = nothing)
    level_code = level_code === nothing ? LevelCode() : level_code
    lead       = ProtoSyn.get_lead(level_code)
    inner_lead = ProtoSyn.get_inner_lead(level_code)

    println(io, lead*"✉  Callback:")
    println(io, inner_lead*"+"*repeat("-", 70)*"+")
    @printf(io, "%s| %-5s | %-25s | %-32s |\n", inner_lead, "Index", "Field", "Value")
    println(io, inner_lead*"+"*repeat("-", 70)*"+")
    @printf(io, "%s| %-5d | %-25s | %-30s   |\n", inner_lead, 1, "Event", cb.event)
    @printf(io, "%s| %-5d | %-25s | %-30s   |\n", inner_lead, 2, "Frequency", cb.frequency)
    println(io, inner_lead*"+"*repeat("-", 70)*"+")
end