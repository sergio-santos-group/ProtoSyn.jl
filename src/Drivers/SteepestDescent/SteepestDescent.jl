module SteepestDescent

using ..Aux
using ..Common
using ..Drivers
using Printf
using LinearAlgebra

@doc raw"""
    Driver(evaluator!::Function[, n_steps::Int64 = 0, f_tol::Float64 = 1e-3, max_step:Float64 = 0.1, callbacks::Tuple{Common.CallbackObject}...])

Define the runtime parameters for the Steepest Descent simulation.
If `n_steps` is zero, a `single point` energy calculation is performed.

# Arguments
- `evaluator!::Function`: Responsible for evaluating the current `state.energy` and calculate the resulting forces. This function should have the following signature:
```
evaluator!(state::Common.State, do_forces::Bool)
```
- `n_steps::Int64`: (Optional) Total amount of steps to be performed (if convergence is not achieved before) (Default: 0).
- `f_tol::Float64`: (Optional) Force tolerance. Defines a finalization criteria, as the steepest descent is considered converged if the maximum force calculated is below this value (Default = 1e-3).
- `max_step::Float64`: (Optional) Defines the maximum value ɣ that the system can jump when applying the forces (Default: 0.1).

# Examples
```julia-repl
julia> Drivers.SteepestDescent.Driver(my_evaluator!)
SteepestDescentDriver(evaluator=my_evaluator!, n_steps=100, f_tol=1e-3, max_step=0.1)

julia> Drivers.SteepestDescent.Driver(my_evaluator!, 0, 1e-6, 0.1)
SteepestDescentDriver(evaluator=my_evaluator!, n_steps=0, f_tol=1e-6, max_step=0.1)
```
!!! tip
    The `my_evaluator!` function often contains an aggregation of pre-defined functions avaliable in [Forcefield](@ref Forcefield). It is possible to combine such functions using the [`@faggregator`](@ref Common) macro.

See also: [`Amber.evaluate!`](@ref Forcefield) [`run!`](@ref)
"""
Base.@kwdef mutable struct DriverConfig{F <: Function}
    
    evaluator!::F
    n_steps::Int64    = 0
    f_tol::Float64    = 1e-3
    max_step::Float64 = 0.1
end
DriverConfig(evaluator!::Function; n_steps::Int64 = 0, f_tol::Float64 = 1e-3, max_step::Float64 = 0.1) = DriverConfig(evaluator!, n_steps, f_tol, max_step)
Base.show(io::IO, b::DriverConfig) = print(io, "SteepestDescent.DriverConfig(evaluator=$(string(b.evaluator!)), n_steps=$(b.n_steps), f_tol=$(b.f_tol), max_step=$(b.max_step))")

#TO DO: Documentation
Base.@kwdef mutable struct DriverState
    
    step::Int64        = 0
    step_size::Float64 = 0.0
    max_force::Float64 = 0.0
    gamma::Float64     = 0.0
end
Base.show(io::IO, b::DriverState) = print(io, "SteepestDescent.DriverState(step=$(b.step), step_size=$(b.step_size), max_force=$(b.max_force), gamma=$(b.gamma))")

# ----------------------------------------------------------------------------------------------------------
#                                                   RUN

@doc raw"""
    run!(state::Common.State, driver::SteepestDescentDriver[, callbacks::::Tuple{Common.CallbackObject}...])

Run the main body of the Driver. If `driver.n_steps` is zero, a `single point` energy calculation is performed.

# Arguments
- `state::Common.State`: Current state of the system to be modified.
- `driver::SteepestDescentDriver`: Defines the parameters for the SteepestDescent simulation. See [`Driver`](@ref).
- `callbacks::Vararg{Common.CallbackObject, N}`: (Optional) Tuple of [`CallbackObject`](@ref Common)s.

The [`CallbackObject`](@ref Common) in this Driver returns the following extra Varargs (in order):
- `max_force::Float64`: The maximum force experienced by the system in the current step.
- `gamma::Float64`: The simulation gamma (ɣ) value for the current step.

!!! tip
    The callback function often contains a [Print](@ref) function.

# Examples
```julia-repl
julia> Drivers.SteepestDescent.run(state, steepest_descent_driver, callback1, callback2, callback3)
```
"""
function run!(state::Common.State, driver_config::DriverConfig, callbacks::Common.CallbackObject...)

    @inline function get_max_force(n_atoms::Int64, forces::Array{Float64, 2})::Float64
        max_force = 0.0
        for i=1:n_atoms
            forceSq = forces[i, 1] ^ 2 + forces[i, 2] ^ 2 + forces[i, 3] ^ 2
            if forceSq > max_force
                max_force = forceSq
            end
        end
        return sqrt(max_force)
    end

    # Evaluate initial energy and forces
    energy::Float64 = driver_config.evaluator!(state, true)
    driver_state = DriverState(0, 0.0, get_max_force(state.size, state.forces), driver_config.max_step)

    # Verify convergence
    if driver_state.max_force < driver_config.f_tol || driver_state.gamma < eps()
        println("⤷ Achieved convergence in $(driver_state.step) steps...")
        return
    end
        
    # Initial callback
    @Common.cbcall callbacks state driver_config driver_state
    
    backup_state = Common.State(state.size)
    while driver_state.step < driver_config.n_steps
        copy!(backup_state, state)
        driver_state.gamma = min(driver_state.gamma, driver_config.max_step)
        driver_state.step_size = driver_state.gamma / get_max_force(state.size, state.forces)
        @. state.xyz += driver_state.step_size * state.forces

        # Calculate new energy and forces
        fill!(state.forces, 0.0)
        energy = driver_config.evaluator!(state, true)
        driver_state.max_force = get_max_force(state.size, state.forces)

        # Verify convergence
        if driver_state.max_force < driver_config.f_tol || driver_state.gamma < eps()
            println("⤷ Achieved convergence in $(driver_state.step) steps...")
            return
        end
        
        # Update gamma
        if energy >= backup_state.energy.eTotal
            driver_state.gamma *= 0.50
            copy!(state, backup_state)
        else
            driver_state.gamma *= 1.05
        end

        driver_state.step += 1
        @Common.cbcall callbacks state driver_config driver_state
    end
end

end