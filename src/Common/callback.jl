export default_energy_step_callback

@doc """
    default_energy_step_callback(n::Int)::Callback

Returns a simple [`Callback`](@ref ProtoSyn.Drivers.Callback) that prints the
current step and total energy value of the [`Pose`](@ref) every `n` steps.

# Examples
```jldoctest
julia> ProtoSyn.Common.default_energy_step_callback(100)
✉  Callback:
+----------------------------------------------------------------------+
| Index | Field                     | Value                            |
+----------------------------------------------------------------------+
| 1     | Event                     | energy_step                      |
| 2     | Frequency                 | 100                              |
+----------------------------------------------------------------------+
```
"""
function default_energy_step_callback(n::Int)::Callback
    function energy_step(pose::Pose, driver_state::DriverState)
        @printf("%8d%12.3f\n", driver_state.step, pose.state.e[:Total])
    end

    return Callback(energy_step, n)
end


export default_energy_step_frame_callback

@doc """
    default_energy_step_frame_callback(n::Int, filename::String)::Callback

Returns a simple [`Callback`](@ref ProtoSyn.Drivers.Callback) that prints the
current step and total energy value of the [`Pose`](@ref) every `n` steps, while
printing the current structure to the given `filename` (in PDB format).

# Examples
```jldoctest
julia> ProtoSyn.Common.default_energy_step_frame_callback(100, "test.pdb")
✉  Callback:
+----------------------------------------------------------------------+
| Index | Field                     | Value                            |
+----------------------------------------------------------------------+
| 1     | Event                     | energy_step_frame                |
| 2     | Frequency                 | 100                              |
+----------------------------------------------------------------------+
```
"""
function default_energy_step_frame_callback(n::Int, filename::String)::Callback
    function energy_step_frame(pose::Pose, driver_state::DriverState)
        @printf("%8d%12.3f\n", driver_state.step, pose.state.e[:Total])
        ProtoSyn.append(pose, filename)
    end

    return Callback(energy_step_frame, n)
end