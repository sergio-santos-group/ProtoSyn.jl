@doc """
    The Common module stores functions with higher levels of abstraction and
    that are common practine in ProtoSyn algorithms.
"""
module Common
    using Printf
    using ProtoSyn
    using ProtoSyn.Drivers: DriverState, Callback
    using ProtoSyn.Units: defaultFloat
    using ProtoSyn.Calculators
    using ProtoSyn.Peptides

    include("energy_function.jl")

    # COMMON CALLBACKS ---------------------------------------------------------

    export default_energy_step_callback

    @doc """
        default_energy_step_callback(n::Int)::Callback

    Returns a simple callback that prints the current step and total energy
    value of the pose every `n` steps.

    # Examples
    ```jldoctest
    julia> default_energy_step_callback(100)
    ProtoSyn.Drivers.Callback(ProtoSyn.Common.var"#energy_step#4"(), 100)
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

    Returns a simple callback that prints the current step and total energy
    value of the pose every `n` steps, while printing the current structures to
    the given `filename` (in PDB format).

    # Examples
    ```jldoctest
    julia> default_energy_step_frame_callback(100)
    ProtoSyn.Drivers.Callback(ProtoSyn.Common.var"#energy_step_frame#4"(), 100)
    ```
    """
    function default_energy_step_frame_callback(n::Int, filename::String)::Callback
        function energy_step_frame(pose::Pose, driver_state::DriverState)
            @printf("%8d%12.3f\n", driver_state.step, pose.state.e[:Total])
            ProtoSyn.append(pose, filename)
        end

        return Callback(energy_step_frame, n)
    end
end