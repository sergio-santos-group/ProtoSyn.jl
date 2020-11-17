mutable struct Callback
    event::Function
    frequency::Int
end

function (callback::Callback)(pose::Pose, driver_state::ProtoSyn.Drivers.DriverState)
    if driver_state.step % callback.frequency == 0
        callback.event(pose, driver_state)
    end
end