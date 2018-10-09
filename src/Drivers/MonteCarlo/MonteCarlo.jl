module MonteCarlo

using ..Common
using Printf

@doc raw"""
    ConfigParameters(n_steps::Int64 = 1000, temperature::Float64 = 0.5)

Define the runtime parameters for the Monte Carlo Driver.

# Arguments
- `n_steps`: Total amount of steps to be performed (Default: 1000).
- `temperature`: Temperature value to be used in the Metropolis Algorithm (Default: 0.5).

# Examples
```julia-repl
julia> Drivers.MonteCarlo.ConfigParameters(1000, 0.5)
Drivers.MonteCarlo.ConfigParameters(n_steps=1000, temperature=0.5)

julia> Drivers.MonteCarlo.ConfigParameters(temperature = 1.2)
Drivers.MonteCarlo.ConfigParameters(n_steps=1000, temperature=1.2)
```
See also: [`load_parameters`](@ref)
"""
struct ConfigParameters

    n_steps::Int64        # (1000)
    temperature::Float64 # (0.5)

    ConfigParameters(; n_steps::Int64 = 1000, temperature::Float64 = 0.5) = new(n_steps, temperature)
end
Base.show(io::IO, b::ConfigParameters) = print(io, "Drivers.MonteCarlo.ConfigParameters(n_steps=$(b.n_steps), temperature=$(b.temperature))")

# ------------------------------------------------------------------------------------------------------------

@doc raw"""
    load_parameters(p::Dict{String, Any})::ConfigParameters

Load the [`ConfigParameters`](@ref) from a dictionary.

# Examples
```julia-repl
julia> Drivers.MonteCarlo.load_parameters(p)
Drivers.MonteCarlo.ConfigParameters(n_steps=1000, temperature=0.5)
```
See also: [`Aux.read_JSON`](@ref Aux)
"""
function load_parameters(p::Dict{String, Any})::ConfigParameters

    return ConfigParameters(p["n_steps"], p["temperature"])

end

# ------------------------------------------------------------------------------------------------------------

@doc raw"""
    run!(state::Common.State, sampler!::Function, evaluator!::Function, params::ConfigParameters[, ostream::IO = stdout, callback::Union{Function, Nothing} = nothing])

Run the main body of the Driver.

# Arguments
- `state::Common.State`: Current state of the system to be modified.
- `sampler!::Function`: Responsible for mutating the current `state`. This function should have the following signature:
```   
sampler!(state::Common.State)
```
- `evaluator!::Function`: Responsible for evaluating the current `state.energy` for the Metropolis Algorithm. This function should have the following signature:
```
evaluator!(state::Common.State, do_forces::Bool)
```
- `params::ConfigParameters`: Hold the runtime parameters of the Monte Carlo Driver. 
- `ostream::IO`: (Optional) Any logging will be written to the supplied `ostream`.
- `callback::Union{Function, Nothing}`: (Optional) If present, this function will be called if the new conformation was accepted. This function should have the following signature:
```
callback(state::Common.State, step::Int64)
```
!!! tip
    The callback function is often a [Print](@ref) function.

# Examples
```julia-repl
julia> Drivers.MonteCarlo.run(state, my_sampler!, my_evaluator!, params!, callback = my_callback)
```
See also: [Mutators](@ref) [`Forcefield.evalenergy!`](@ref Forcefield)
"""
function run!(
    state::Common.State,
    sampler!::Function,
    evaluator!::Function,
    params::ConfigParameters;
    ostream::IO = stdout,
    callback::Union{Function, Nothing} = nothing)

    xyz0 = copy(state.xyz)
    ene0 = evaluator!(state, false)

    step::Int64 = 0
    while step < params.n_steps
        step += 1

        #Generate new state and evaluate its energy
        sampler!(state)
        ene1 = evaluator!(state, false)

        if ene1 < ene0 || rand() < exp(-(ene1 - ene0) / params.temperature)
            ene0 = ene1
            xyz0[:] = state.xyz
            if callback != nothing
                callback(state, step)
            end
            write(ostream, @sprintf "(MC) Step: %4d | Energy: %9.4f\n" step ene1)
        else
            state.xyz[:] = xyz0
        end
    end
end

end