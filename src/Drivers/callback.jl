mutable struct Callback
    event::Function # Receives a pose and a DriverState
    frequency::Int
end

function (callback::Callback)(pose::Pose, driver_state::ProtoSyn.Drivers.DriverState)
    if driver_state.step % callback.frequency == 0
        callback.event(pose, driver_state)
    end
end