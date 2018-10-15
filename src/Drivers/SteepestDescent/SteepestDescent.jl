module SteepestDescent

using ..Common
using Printf
using LinearAlgebra

@doc raw"""
    ConfigParameters(n_steps::Int64 = 0, log_freq::Int64 = 0, f_tol::Float64 = 1e-3, max_step:Float64 = 0.1)

Define the runtime parameters for the Monte Carlo Driver.
If `n_steps` is zero, a `single point` energy calculation is performed.

# Arguments
- `n_steps`: Total amount of steps to be performed (if convergence is not achieved before) (Default: 0).
- `log_freq`: Defines the frequency (in steps) of the logs (Default: 0).
- `f_tol`: Force tolerance. Defines a finalization criteria, as the steepest descent is considered converged if the maximum force calculated is below this value (Default = 1e-3).
- `max_step`: Defines the maximum value ɣ that the system can jump when applying the forces (Default: 0.1).

# Examples
```julia-repl
julia> Drivers.SteepestDescent.ConfigParameters(100, 5, 1e-3, 0.1)
Drivers.SteepestDescent.ConfigParameters(n_steps=100, log_freq=5, f_tol=1e-3, max_step=0.1)

julia> Drivers.SteepestDescent.ConfigParameters(f_tol = 1e-6)
Drivers.SteepestDescent.ConfigParameters(n_steps=0, log_freq=0, f_tol=1e-6, max_step=0.1)
```
See also: [`load_parameters`](@ref)
"""
mutable struct ConfigParameters
    
    n_steps::Int64
    log_freq::Int64
    callback_freq::Int64
    f_tol::Float64
    max_step::Float64

    ConfigParameters(; n_steps::Int64 = 0, log_freq::Int64 = 1, callback_freq::Int64 = 0, f_tol::Float64 = 1e-3, max_step::Float64 = 0.1) = new(n_steps, log_freq, callback_freq, f_tol, max_step)
end
Base.show(io::IO, b::ConfigParameters) = print(io, "Drivers.SteepestDescent.ConfigParameters(n_steps=$(b.n_steps), log_freq=$(b.log_freq), f_tol=$(b.f_tol), max_step=$(b.max_step))")

# ------------------------------------------------------------------------------------------------------------

@doc raw"""
    load_parameters(p::Dict{String, Any})::ConfigParameters

Load the [`ConfigParameters`](@ref) from a dictionary.

# Examples
```julia-repl
julia> Drivers.SteepestDescent.load_parameters(p)
Drivers.SteepestDescent.ConfigParameters(n_steps=100, log_freq=5, f_tol=1e-3, max_step=0.1)
```
See also: [`Aux.read_JSON`](@ref Aux)
"""
function load_parameters(p::Dict{String, Any})::ConfigParameters

    return ConfigParameters(p["n_steps"], p["log_freq"], p["f_tol"], p["max_step"])

end

#TO DO:
# GAMMA SHOULD BE UPDATED BASED ON ENERGY (?)
@doc raw"""
    run!(state::Common.State, evaluator!::Function, params::ConfigParameters[, ostream::IO = stdout, callback::Union{Function, Nothing} = nothing])

Run the main body of the Driver. If `params.n_steps` is zero, a `single point` energy calculation is performed.

# Arguments
- `state::Common.State`: Current state of the system to be modified.
- `evaluator!::Function`: Responsible for evaluating the current `state.energy` and calculate the resulting forces. This function should have the following signature:
```
evaluator!(state::Common.State, do_forces::Bool)
```
- `params::ConfigParameters`: Hold the runtime parameters of the Steepest Descent Driver. 
- `ostream::IO`: (Optional) Any logging will be written to the supplied `ostream`.
- `callback::Union{Function, Nothing}`: (Optional) If present, this function will be called every [`ConfigParameters`](@ref).log_freq. This function should have the following signature:
```
callback(state::Common.State, step::Int64)
```
!!! tip
    The callback function is often a [Print](@ref) function.

# Examples
```julia-repl
julia> Drivers.SteepestDescent.run(state, my_evaluator!, params!, callback = my_callback)
```
See also: [`Forcefield.evalenergy!`](@ref Forcefield)
"""
function run!(
    state::Common.State,
    evaluator!::Function,
    params::ConfigParameters;
    ostream::IO = stdout,
    callback::Union{Function, Nothing} = nothing)

    @inline function get_max_force(f::Array{Float64, 2})
        return sqrt(maximum(sum(f.*f, dims = 2)))
    end

    @inline function logger(s::Int64, e::Float64, f::Float64, g::Float64)
        write(ostream, @sprintf "(SD) Step: %4d | Energy: %9.4f | maxForce: %9.4f | gamma: %8.3e\n" s e f g)
    end
    
    @inline function is_output_step(freq::Int64, step::Int64)
        return ((freq>0) && (step%freq==0))
    end


    # Evaluate initial energy and forces
    energy::Float64 = evaluator!(state, true)
    max_force::Float64  = get_max_force(state.forces)
    energy_old = energy
    
    gamma::Float64 = 1.0#e-5
    step::Int64 = 0
    
    
    if callback != nothing
        callback(state, step)
    end
    logger(step, energy_old, max_force, gamma)
    
    while step < params.n_steps
        step += 1

        #Update system coordinates
        gamma = min(gamma, params.max_step)
        stepsize = gamma / get_max_force(state.forces)
        @. state.xyz += stepsize * state.forces

        # housekeep variables
        energy_old = energy

        # Calculate new energy and forces
        fill!(state.forces, 0.0)
        energy = evaluator!(state, true)
        max_force = get_max_force(state.forces)

        # call callback function and output information to log 
        if is_output_step(params.callback_freq, step) && (callback != nothing)
            callback(state, step)
        end

        if is_output_step(params.log_freq, step)
            logger(step, energy, max_force, gamma)
        end
        
        # check if force convergence was achieved
        if max_force < params.f_tol
            write(ostream, "Achieved convergence (f_tol < $(params.f_tol)) in $step steps.\n")
            break
        end
        
        # check if gamma is below machine precision
        if gamma < eps()
            write(ostream, "Gamma below machine precision! Exiting ...\n")
            break
        end
        
        # Update gamma
        if energy > energy_old
            gamma *= 0.50
        else
            gamma *= 1.05
        end
    end


    if callback != nothing
        callback(state, step)
    end
    logger(step, energy, max_force, gamma)
end

end