module SteepestDescent

using ..Aux
using ..Common
using ..Drivers

@doc raw"""
    Driver(evaluator!::Function[, n_steps::Int64 = 0, f_tol::Float64 = 1e-3, max_displacement:Float64 = 0.1, callbacks::Tuple{Common.CallbackObject}...])

Define the runtime parameters for the Steepest Descent simulation.
If `n_steps` is zero, a `single point` energy calculation is performed.

# Arguments
- `evaluator!::Function`: Responsible for evaluating the current `state.energy` and calculate the resulting forces. This function should have the following signature:
```
evaluator!(state::Common.State, do_forces::Bool)
```
- `n_steps::Int64`: (Optional) Total amount of steps to be performed (if convergence is not achieved before) (Default: 0).
- `f_tol::Float64`: (Optional) Force tolerance. Defines a finalization criteria, as the steepest descent is considered converged if the maximum force calculated is below this value (Default = 1e-3).
- `max_displacement::Float64`: (Optional) Defines the maximum value ɣ that the system can jump when applying the forces (Default: 0.1).

# Examples
```julia-repl
julia> Drivers.SteepestDescent.Driver(my_evaluator!)
SteepestDescentDriver(evaluator=my_evaluator!, n_steps=100, f_tol=1e-3, max_displacement=0.1)

julia> Drivers.SteepestDescent.Driver(my_evaluator!, 0, 1e-6, 0.1)
SteepestDescentDriver(evaluator=my_evaluator!, n_steps=0, f_tol=1e-6, max_displacement=0.1)
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

DriverConfig(evaluator!::F) where {F <: Function} = DriverConfig(evaluator! = evaluator!)

function Base.show(io::IO, b::DriverConfig)
    print(io, "SteepestDescent.DriverConfig")
    for p in fieldnames(DriverState)
        print(io, "\n   $(String(p)) = $(getproperty(b,p))")
    end
end


#TO DO: Documentation
Base.@kwdef mutable struct DriverState <: Drivers.AbstractDriverState
    
    step::Int64        = 0
    step_size::Float64 = -1.0
    max_force::Float64 = -1.0
    converged::Bool    = false
    stalled::Bool      = false
end

function Base.show(io::IO, b::DriverState)
    print(io, "SteepestDescent.DriverState")
    for p in fieldnames(DriverState)
        print(io, "\n   $(String(p)) = $(getproperty(b,p))")
    end
end


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

    # Evaluate initial energy and forces
    #   start by calculating nonbonded lists:
    #   assuming that each atomic displacement is limited
    #   to 'max_displacement' at each step, then any
    #   2 atoms, if moving in oposite directions, can be
    #   displaced 2*max_displacement from each other.
    #   Therefore, the maximum total displacement during
    #   'nblist_freq' steps is 2*max_displacement*nblist_freq.
    #   Hence, the nblist cutoff is (cutoff + 2*max_displacement*nblist_freq).
    if state.nblist != nothing && state.nblist.cutoff > 0.0
        state.nblist.buffer  = 2.0 * driver_config.max_displacement
        state.nblist.buffer *= driver_config.nblist_freq
    end
    
    update_nblist(state)
    energy = driver_config.evaluator!(state, true)

    # instantiate a new DriverState object.
    # By default, no optimization step has yet been taken
    # apart from calculating the energy and forces for the
    # input state
    driver_state = DriverState()
    driver_state.max_force = get_max_force(state.size, state.forces)
    
    # if max force is already below the requested force tolerance,
    # simply return the current driver state
    driver_state.converged = driver_state.max_force < driver_config.f_tol

    if driver_state.converged
        return driver_state
    end
    
    # Initial callbacks
    @Common.cbcall callbacks state driver_state driver_config
    
    # Initialize γ and backup copy
    γ = 1.0
    backup_state = Common.State(state.size)
    
    #region MAINLOOP
    while driver_state.step < driver_config.n_steps
        
        # save current state to backup
        copy!(backup_state, state)

        # calculate current scaling factor and step size
        γ = min(γ, driver_config.max_displacement)
        driver_state.step_size = γ / driver_state.max_force

        # update coordinates
        @. state.xyz += driver_state.step_size * state.forces

        # update nonbonded lists if required
        # TODO: check nblists 
        if (driver_config.nblist_freq > 0) &&
            (driver_state.step > 0) &&
            (driver_state.step % driver_config.nblist_freq == 0)
            update_nblist(state)
        end

        # Calculate new energy and forces
        # (make sure to reset forces)
        fill!(state.forces, 0.0)
        energy = driver_config.evaluator!(state, true)
        driver_state.max_force = get_max_force(state.size, state.forces)

        # Verify convergence
        driver_state.converged = driver_state.max_force < driver_config.f_tol
        driver_state.stalled = γ < eps()

        if driver_state.converged || driver_state.stalled
            break
        end
        
        # Update gamma
        if energy >= backup_state.energy.total
            γ *= 0.50
            copy!(state, backup_state)
            driver_state.max_force = get_max_force(state.size, state.forces)
        else
            γ *= 1.05
        end

        # update current step and call calback functions (if any)
        driver_state.step += 1
        @Common.cbcall callbacks state driver_state driver_config

    end
    #endregion

    return driver_state
end

end