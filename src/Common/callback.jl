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


"""
# TODO DOCUMENTATION
"""
function default_energy_step_detailed(n::Int, msg::String = "Callback", color::Symbol = :none)
    function energy_step_detailed(pose::Pose, driver_state::ProtoSyn.Drivers.DriverState)
        N = 20 + (20 * length(pose.state.e)) - length(msg)
        M = floor(Int, N/2) - 3
        E = N-6-(M*2)
        acc = driver_state.step === 0 ? 1.0 : driver_state.acceptance_count/driver_state.step
        if driver_state.step === 0
            s  = @sprintf("%s\n", "| "*repeat("-", M)*" "*repeat(" ", E)*msg*" "*repeat("-", M)*" |")
            s *= @sprintf("%6s%14s", "Step", "Accept-Ratio")
            for component in sort(collect(keys(pose.state.e)))
                s *= @sprintf("%20s", string(component))
            end
            s*= "\n"
            printstyled(s, color = :cyan)
        end

        s = @sprintf("%6d%14.3f", driver_state.step, acc)
        for component in sort(collect(keys(pose.state.e)))
            s *= @sprintf("%20.3f", pose.state.e[component])
        end
        s*= "\n"
        printstyled(s, color = color)
    end

    return Callback(energy_step_detailed, n)
end