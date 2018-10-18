module SteepestDescent

using ..Aux
using ..Common
using Printf
using LinearAlgebra

@doc raw"""
    SteepestDescentDriver(evaluator!::Function[, n_steps::Int64 = 0, f_tol::Float64 = 1e-3, max_step:Float64 = 0.1])

Define the runtime parameters for the Monte Carlo simulation.
If `n_steps` is zero, a `single point` energy calculation is performed.

# Arguments
- `evaluator!::Function`: Responsible for evaluating the current `state.energy` and calculate the resulting forces. This function should have the following signature:
```
evaluator!(state::Common.State, do_forces::Bool)
```
- `n_steps`: (Optional) Total amount of steps to be performed (if convergence is not achieved before) (Default: 0).
- `f_tol`: (Optional) Force tolerance. Defines a finalization criteria, as the steepest descent is considered converged if the maximum force calculated is below this value (Default = 1e-3).
- `max_step`: (Optional) Defines the maximum value ɣ that the system can jump when applying the forces (Default: 0.1).
- `ostream`: (Optional) Defines the output stream for logging

# Examples
```julia-repl
julia> Drivers.SteepestDescent.ConfigParameters(Forcefield.evalenergy!, 100, 1e-3, 0.1)
Drivers.SteepestDescent.ConfigParameters(evaluator!=Forcefield.evalenergy!, n_steps=100, f_tol=1e-3, max_step=0.1)

julia> Drivers.SteepestDescent.ConfigParameters(Forcefield.evalenergy!, f_tol = 1e-6)
Drivers.SteepestDescent.ConfigParameters(evaluator!=Forcefield.evalenergy!, n_steps=0, f_tol=1e-6, max_step=0.1)
```
See also: [`load_parameters`](@ref) [`Forcefield.Amber.evalenergy!`](@ref Forcefield)
"""
mutable struct SteepestDescentDriver
    
    evaluator!::Function
    n_steps::Int64
    f_tol::Float64
    max_step::Float64

    SteepestDescentDriver(evaluator!::Function; n_steps::Int64 = 0, f_tol::Float64 = 1e-3, max_step::Float64 = 0.1) = new(evaluator!, n_steps, f_tol, max_step)
end
Base.show(io::IO, b::SteepestDescentDriver) = print(io, "SteepestDescentDriver(evaluator!=$(string(b.evaluator!)), n_steps=$(b.n_steps), f_tol=$(b.f_tol))")

# ----------------------------------------------------------------------------------------------------------
#                                                   RUN

@doc raw"""
    run!(state::Common.State, driver::SteepestDescentDriver[, callback::Union{Common.CallbackObject, Nothing} = nothing])

Run the main body of the Driver. If `driver.n_steps` is zero, a `single point` energy calculation is performed.

# Arguments
- `state::Common.State`: Current state of the system to be modified.
- `driver::SteepestDescentDriver`: Defines the parameters for the SteepestDescent simulation. See [`SteepestDescentDriver`](@ref).
- `callbacks::Vararg{Common.CallbackObject, N}`: (Optional) Tuple of [`CallbackObject`](@ref)s (Default: empty).

!!! tip
    The callback function often contains a [Print](@ref) function.

# Examples
```julia-repl
julia> Drivers.SteepestDescent.run(state, steepest_descent_driver, callback1, callback2, callback3)
```
"""
function run!(state::Common.State, driver::SteepestDescentDriver, callbacks::Common.CallbackObject...)

    println("CALLBACKS:", callbacks)

    @inline function get_max_force(f::Array{Float64, 2})
        return sqrt(maximum(sum(f.*f, dims = 2)))
    end

    # Evaluate initial energy and forces
    energy::Float64 = evaluator!(state, true)
    max_force::Float64  = get_max_force(state.forces)
    energy_old = energy
    
    gamma::Float64 = 1.0
    step::Int64 = 0
    
    # Initial callback
    # @Common.cbcall callbacks step state driver
    
    while step < params.n_steps
        step += 1

        # Update system coordinates
        gamma = min(gamma, params.max_step)
        stepsize = gamma / get_max_force(state.forces)
        @. state.xyz += stepsize * state.forces

        # Housekeep variables
        energy_old = energy

        # Calculate new energy and forces
        fill!(state.forces, 0.0)
        energy = evaluator!(state, true)
        max_force = get_max_force(state.forces)

        # Call callback function and output information to log 
        # @Common.cbcall callbacks step state driver
        
        # Check if force convergence was achieved
        if max_force < params.f_tol
            write(driver.ostream, "Achieved convergence (f_tol < $(driver.f_tol)) in $step steps.\n")
            break
        end
        
        # Check if gamma is below machine precision
        if gamma < eps()
            write(driver.ostream, "Gamma below machine precision! Exiting ...\n")
            break
        end
        
        # Update gamma
        if energy > energy_old
            gamma *= 0.50
        else
            gamma *= 1.05
        end
    end

    #Final callback
    # @Common.cbcall callbacks step state driver
end

end