module SteepestDescent

using ..Aux
using ..Common
using ..Drivers
using ..Abstract

#= ----------------------------------------------------------------------------
Steepest Descent (or Gradient Descent)

# Algorithm explanation:
First-order iterative optimization algorithm, calculated the forces applying
in the system as the first derivative of the energy function, and applies
the calculated forces acting over each atom in a step wise fashion.
The `step_size` is updated each cycle based on the `max_force` being felt by
the system (smaller the bigger the force), up until a defined `max_step`.

# Convergence criteria
1. Number of `n_steps` performed;
2. Max force is below the defined tolerance `f_tol`;
3. Gamma (γ) is below the machine precision (stalled);
---------------------------------------------------------------------------- =#

@doc raw"""
    DriverConfig(; evaluator::Abstract.Evaluator, n_steps::Int64 = 0, nblist_freq:::Int64 = 0, f_tol::Float64 = 1e-3, max_step:Float64 = 0.1, callbacks::Vector{Abstract.CallbackObject} = []])

Define the runtime parameters for the Steepest Descent simulation.
If `n_steps` is zero, a `single point` energy calculation is performed.

# Arguments
- `evaluator!::Abstract.Evaluator`: Responsible for evaluating the current `state.energy` and calculate the resulting forces. See [`Evaluator`](@ref Forcefield)
- `nblist_freq::Int64`: (Optional) Frequency to update neighbouring lists in the current [`State`](@ref Common) (Default: 0).
- `n_steps::Int64`: (Optional) Total amount of steps to be performed (if convergence is not achieved before) (Default: 0).
- `f_tol::Float64`: (Optional) Force tolerance. Defines a finalization criteria, as the steepest descent is considered converged if the maximum force calculated is below this value (Default = 1e-3).
- `max_step::Float64`: (Optional) Defines the maximum value ɣ that the system can jump when applying the forces (Default: 0.1).
- `callbacks::Vector{Abstract.CallbackObject}` (Optional) List of [`CallbackObject`](@ref Common)s (Default: empty)

# Examples
```julia-repl
julia> Drivers.SteepestDescent.Driver(evaluator = my_evaluator!)
SteepestDescent.DriverConfig
   evaluator = (...)
   n_steps = 0
   nblist_freq = 0
   f_tol = 0.001
   max_step = 0.1
   callbacks = ProtoSyn.Common.CallbackObject[]
```
!!! tip
    The `my_evaluator!` function often contains an aggregation of pre-defined functions avaliable in [Forcefield](@ref Forcefield).

See also: [`Amber.evaluate!`](@ref Forcefield) [`run!`](@ref)
"""
mutable struct DriverConfig <: Abstract.DriverConfig
    
    evaluator::Abstract.Evaluator            # Required
    n_steps::Int64                           # Default: 0
    nblist_freq::Int64                       # Default: 0
    f_tol::Float64                           # Default: 1e-3
    max_step::Float64                        # Default: 0.1
    callbacks::Vector{Common.CallbackObject} # Default: empty

    DriverConfig(; evaluator::Abstract.Evaluator,
        n_steps::Int64 = 0,
        nblist_freq::Int64 = 0,
        f_tol::Float64 = 1e-3,
        max_step::Float64 = 0.1,
        callbacks::Vector{<:Abstract.CallbackObject} = Vector{Common.CallbackObject}()) = begin
            new(evaluator, n_steps, nblist_freq, f_tol, max_step, callbacks)
    end
end

function Base.show(io::IO, b::DriverConfig)
    print(io, "SteepestDescent.DriverConfig")
    for p in fieldnames(DriverConfig)
        print(io, "\n   $(String(p)) = $(getproperty(b,p))")
    end
end


@doc raw"""
    DriverState(;step::Int64 = 0, step_size::Float64 = -1.0, max_force::Float64 = -1.0, converged::Bool = false, stalled::Bool = false)

Define the runtime parameters given by the Steepest Descent simulation.

# Arguments
- `step::Int64`: Current step of the simulation.
- `step_size::Float64`: Current step_size of the simulation.
- `max_force::Float64`: Current max_force of the system.
- `converged::Bool`: A system is considered converged in the `max_force` value is below the defined `f_tol` value in [`DriverConfig`](@ref)
- `stalled::Bool`: A system is considered stalled if the gamma (γ) value is below machine precision.

# Examples
```julia-repl
julia> Drivers.SteepestDescent.DriverState()
SteepestDescent.DriverState
   step = 0
   step_size = -1.0
   max_force = -1.0
   converged = false
   stalled = false

```
See also: [`run!`](@ref)
"""
Base.@kwdef mutable struct DriverState <: Abstract.DriverState
    
    # Parameter:       # Default:
    step::Int64        = 0
    step_size::Float64 = -1.0
    max_force::Float64 = -1.0
    converged::Bool    = false
    stalled::Bool      = false
end


# ----------------------------------------------------------------------------------------------------------
#                                                   RUN

@doc raw"""
    run!(state::Common.State, driver_config::DriverConfig)

Run the main body of the Driver. If `driver_config.n_steps` is zero, a `single point` energy calculation is performed.

# Arguments
- `state::Common.State`: Current state of the system to be modified.
- `driver_config::DriverConfig`: Defines the parameters for the SteepestDescent simulation. See [`DriverConfig`](@ref).

# Examples
```julia-repl
julia> Drivers.SteepestDescent.run!(state, steepest_descent_driver)
```
"""
function run!(state::Common.State, driver_config::DriverConfig)

    # Evaluate initial energy and forces
    #   start by calculating nonbonded lists:
    #   assuming that each atomic displacement is limited
    #   to 'max_displacement' at each step, then any
    #   2 atoms, if moving in oposite directions, can be
    #   displaced 2*max_displacement from each other.
    #   Therefore, the maximum total displacement during
    #   'nblist_freq' steps is 2*max_displacement*nblist_freq.
    #   Hence, the nblist cut_off is (cut_off + 2*max_displacement*nblist_freq).
    if state.nblist != nothing && state.nblist.cut_off > 0.0
        state.nblist.buffer  = 2.0 * driver_config.max_displacement
        state.nblist.buffer *= driver_config.nblist_freq
    end # end if
    
    Common.update_nblist!(state)
    energy = driver_config.evaluator.evaluate!(state, driver_config.evaluator.components, true)

    # instantiate a new DriverState object.
    # By default, no optimization step has yet been taken
    # apart from calculating the energy and forces for the
    # input state
    driver_state = DriverState()
    driver_state.max_force = Aux.get_max_force(state.size, state.forces)
    
    # if max force is already below the requested force tolerance,
    # simply return the current driver state
    driver_state.converged = driver_state.max_force < driver_config.f_tol

    if driver_state.converged
        return driver_state
    end # end if
    
    # Initial callbacks
    Common.@cbcall driver_config.callbacks state driver_state
    
    # Initialize γ and backup copy
    γ = 1.0
    backup_state = Common.State(state.size)
    
    #region MAINLOOP
    while driver_state.step < driver_config.n_steps
        
        # save current state to backup
        copy!(backup_state, state)

        # calculate current scaling factor and step size
        γ = min(γ, driver_config.max_step)
        driver_state.step_size = γ / driver_state.max_force

        # update coordinates
        @. state.xyz += driver_state.step_size * state.forces

        # update nonbonded lists if required
        # TODO: check nblists 
        if (driver_config.nblist_freq > 0) &&
            (driver_state.step > 0) &&
            (driver_state.step % driver_config.nblist_freq == 0)
            Common.update_nblist!(state)
        end # end if

        # Calculate new energy and forces
        # (make sure to reset forces)
        fill!(state.forces, 0.0)
        energy = driver_config.evaluator.evaluate!(state, driver_config.evaluator.components, true)
        driver_state.max_force = Aux.get_max_force(state.size, state.forces)

        # Verify convergence
        driver_state.converged = driver_state.max_force < driver_config.f_tol
        driver_state.stalled = γ < eps()

        if driver_state.converged || driver_state.stalled
            break
        end # end if
        
        # Update gamma
        if energy >= backup_state.energy.total
            γ *= 0.50
            copy!(state, backup_state)
            driver_state.max_force = Aux.get_max_force(state.size, state.forces)
        else
            γ *= 1.05
        end # end if

        # update current step and call calback functions (if any)
        driver_state.step += 1
        @Common.cbcall driver_config.callbacks state driver_state
    end # end while
    #endregion

    return driver_state
end # end function

end # end module