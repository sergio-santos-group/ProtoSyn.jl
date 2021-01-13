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

    # COMMON ENERGY FUNCTIONS --------------------------------------------------

    export default_energy_function

    @doc """
        default_energy_function(::Type{T}) where {T <: AbstractFloat}
        default_energy_function()

    Returns a default energy function for ProtoSyn. As of 2021, this includes
    the following terms:
    +----------------------------------------------------------+
    | Index | Component name                      | Weight (ɑ) |
    +----------------------------------------------------------+
    | 1     | Caterpillar_Solvation               | 0.01       |
    | 2     | TorchANI_ML_Model                   | 1.00       |
    +----------------------------------------------------------+
    **Note:** If no Type{T} is provided, will use ProtoSyn.Units.defaultFloat;

    # Examples
    ```jldoctest
    julia> ProtoSyn.Common.default_energy_function()
    ...
    ```
    """
    default_energy_function(::Type{T}) where {T <: AbstractFloat} = begin
        return Calculators.EnergyFunction(Dict(
            Calculators.TorchANI.torchani_model => T(1.0),
            Peptides.Calculators.Caterpillar.solvation_energy => T(0.01),
            Calculators.Restraints.bond_distance_restraint => T(0.05)
        ))
    end

    default_energy_function() = begin
        default_energy_function(defaultFloat)
    end

    export default_energy_step_callback

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
        default_energy_step_frame_callback(n::Int)::Callback

    Returns a simple callback that prints the current step and total energy
    value of the pose every `n` steps, while printing the current structures to
    the given filename (in PDB format).

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