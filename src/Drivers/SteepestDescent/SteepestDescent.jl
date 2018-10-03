module SteepestDescent

using ..Common
using Printf

@doc raw"""
    ConfigParameters(n_steps::Int64, log_freq::Int64, f_tol::Float64, max_step:Float64)

Define the runtime parameters for the Monte Carlo Driver.

# Arguments
- `n_steps`: Total amount of steps to be performed (if convergence is not achieved before).
- `log_freq`: Defines the frequency (in steps) of the logs.
- `f_tol`: Force tolerance. Defines a finalization criteria, as the steepest descent is considered converged if the maximum force calculated is below this value.
- `max_step`: Defines the maximum value ɣ that the system can jump when applying the forces.

# Examples
```julia-repl
julia> Drivers.SteepestDescent.ConfigParameters(100, 5, 1e-3, 0.1)
Drivers.SteepestDescent.ConfigParameters(n_steps=100, log_freq=5, f_tol=1e-3, max_step=0.1)
```
See also: [`load_parameters`](@ref)
"""
struct ConfigParameters
    
    n_steps::Int64
    log_freq::Int64
    f_tol::Float64
    max_step::Float64

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

Run the main body of the Driver.

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

    function get_max_force(f::Array{Float64, 2})
        return sqrt(maximum(sum(f.^2, dims = 2)))
    end

    # Evaluate initial energy and forces
    energy_old::Float64 = evaluator!(state, true)
    max_force = get_max_force(state.forces)
    old_max_force::Float64 = max_force
    
    
    gamma::Float64 = 1e-5
    step::Int64 = 0
    if callback != nothing
        callback(state, step)
    end
    while step < params.n_steps
        step += 1

        #Update system coordinates
        gamma = min(gamma, params.max_step)
        dx = (gamma / get_max_force(state.forces)) * state.forces
        state.xyz += dx

        #Calculate new energy and forces
        state.forces = zeros(size(state.xyz, 1), 3)
        energy = evaluator!(state, true)
        max_force = get_max_force(state.forces)

        #Callback, if present
        if params.log_freq != 0 && step % params.log_freq == 0
            if callback != nothing
                callback(state, step)
            end
            write(ostream, @sprintf "(SD) Step: %4d | Energy: %9.4f | maxForce: %9.4f\n" step energy max_force)
        end

        #If max force variance is bellow fTol threshold, convergence was achieved: exit function
        if abs(old_max_force - max_force) < params.f_tol
            if callback != nothing
                callback(state, step)
            end
            write(ostream, "Achieved convergence (below fTol $(params.f_tol)) in $step steps.\n")
            break
        end

        #Update gamma constant
        if max_force > old_max_force
            gamma *= 0.50
        else
            gamma *= 1.05
        end

        #Update current force to new cycle
        old_max_force = max_force
    end
end

end