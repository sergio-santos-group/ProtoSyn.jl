using Printf

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
    println(" Callback:")
    println(io, "+"*repeat("-", 58)*"+")
    @printf(io, "| %-5s | %-15s | %-30s |\n", "Index", "Field", "Value")
    println(io, "+"*repeat("-", 58)*"+")
    @printf(io, "| %-5d | %-15s | %-30d |\n", 1, "Frequency", cb.frequency)
    @printf(io, "| %-5d | %-15s | %-30s |\n", 2, "Event", "$(cb.event)")
    println(io, "+"*repeat("-", 58)*"+")
end