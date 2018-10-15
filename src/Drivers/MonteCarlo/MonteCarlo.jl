module MonteCarlo

using ..Aux
using ..Common
using Printf

#TODO: Document structure
mutable struct MonteCarloDriver
    sampler! :: Function
    evaluator! :: Function
    temperature :: Float64
    nsteps::Int64
end

#TODO: Update documentation
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
function run!(state::Common.State, driver::MonteCarloDriver, callback::Union{Common.CallbackObject, Nothing}=nothing)
    
    step = 0
    xyz0 = copy(state.xyz)
    ene0 = driver.evaluator!(state, false)
    acceptance_count = 0

    while step < driver.nsteps
        step += 1
        driver.sampler!(state)
        ene1 = driver.evaluator!(state, false)

        if (ene1 < ene0) || (rand() < exp(-(ene1 - ene0) / driver.temperature))
            ene0 = ene1
            xyz0[:] = state.xyz
            acceptance_count += 1
        else
            state.xyz[:] = xyz0
        end
        
        @Aux.cbcall callback step state driver (acceptance_count/step)
    end
end

end