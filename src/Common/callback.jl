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
printing the current structure to the given `filename` (format infered from
`filename` extension).

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


@doc """
    default_energy_step_detailed(n::Int, msg::String = "Callback", color::Symbol = :none, file_out::Opt{String} = nothing, print_to_sdtout::Bool = true)::Callback

Returns a detailed [`Callback`](@ref ProtoSyn.Drivers.Callback) that prints the
current step and all energy components of the [`Pose`](@ref) every `n` steps.
An optional `msg` String can be provided as a header. Output can be colored by
setting the `color` flag. For a list of available colors, check the
`printstyled` documentation. By default, `default_energy_step_detailed` prints
the output to the stdout. This can be toggled by setting the `print_to_sdtout`
flag. Optionally, if a `file_out` is provided (as a String), output is also
printed to that file.

# Examples
```jldoctest
julia> ProtoSyn.Common.default_energy_step_detailed(1)
✉  Callback:
+----------------------------------------------------------------------+
| Index | Field                     | Value                            |
+----------------------------------------------------------------------+
| 1     | Event                     | energy_step_detailed             |
| 2     | Frequency                 | 1                                |
+----------------------------------------------------------------------+
```
"""
function default_energy_step_detailed(n::Int, msg::String = "Callback", color::Symbol = :none, file_out::Opt{String} = nothing, print_to_sdtout::Bool = true)::Callback
    function energy_step_detailed(pose::Pose, driver_state::ProtoSyn.Drivers.DriverState)
        if file_out !== nothing
            io = open(file_out, "a")
        end
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
            print_to_sdtout && printstyled(s, color = :cyan)
            if file_out !== nothing
                print(io, s)
            end
        end

        s = @sprintf("%6d%14.3f", driver_state.step, acc)
        for component in sort(collect(keys(pose.state.e)))
            s *= @sprintf("%20.3f", pose.state.e[component])
        end
        s*= "\n"
        print_to_sdtout && printstyled(s, color = color)
        if file_out !== nothing
            print(io, s)
            close(io)
        end
    end

    return Callback(energy_step_detailed, n)
end


@doc """
default_energy_step_frame_detailed(n::Int, output_frame::String, msg::String = "Callback", color::Symbol = :none, output_log::Opt{String} = nothing, print_to_sdtout::Bool = true)::Callback

Returns a detailed [`Callback`](@ref ProtoSyn.Drivers.Callback) that prints the
current step and all energy components of the [`Pose`](@ref) every `n` steps,
while printing the current structure to the given `filename` (format infered
from `filename` extension). An optional `msg` String can be provided as a
header. Output can be colored by setting the `color` flag. For a list of
available colors, check the `printstyled` documentation. By default,
`default_energy_step_detailed` prints the output to the stdout. This can be
toggled by setting the `print_to_sdtout` flag. Optionally, if a `file_out` is
provided (as a String), output is also printed to that file.

# Examples
```jldoctest
julia> ProtoSyn.Common.default_energy_step_detailed(1)
✉  Callback:
+----------------------------------------------------------------------+
| Index | Field                     | Value                            |
+----------------------------------------------------------------------+
| 1     | Event                     | energy_step_detailed             |
| 2     | Frequency                 | 1                                |
+----------------------------------------------------------------------+
```
"""
function default_energy_step_frame_detailed(n::Int, output_frame::String, msg::String = "Callback", color::Symbol = :none, output_log::Opt{String} = nothing, print_to_sdtout::Bool = true)::Callback
    function energy_step_frame_detailed(pose::Pose, driver_state::ProtoSyn.Drivers.DriverState)
        
        # Write structure to output file
        ProtoSyn.append(pose, output_frame)

        # Write log to output file
        if output_log !== nothing
            io = open(output_log, "a")
        end
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
            print_to_sdtout && printstyled(s, color = :cyan)
            if output_log !== nothing
                print(io, s)
            end
        end

        s = @sprintf("%6d%14.3f", driver_state.step, acc)
        for component in sort(collect(keys(pose.state.e)))
            s *= @sprintf("%20.3f", pose.state.e[component])
        end
        s*= "\n"
        print_to_sdtout && printstyled(s, color = color)
        if output_log !== nothing
            print(io, s)
            close(io)
        end
    end

    return Callback(energy_step_frame_detailed, n)
end