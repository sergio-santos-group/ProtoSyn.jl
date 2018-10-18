module MonteCarlo

using ..Aux
using ..Common
using Printf

@doc raw"""
    MonteCarloDriver(sampler!::Function, evaluator!::Function, [, temperature::Float64 = 1.0, n_steps::Int64 = 1])
Define the runtime parameters for the Monte Carlo simulation.
If `n_steps` is one, a `single point` Monte Carlo movement is performed and either accepted or not.

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
mutable struct MonteCarloDriver

    sampler! :: Function
    evaluator! :: Function
    temperature :: Float64
    n_steps::Int64

end
MonteCarloDriver(sampler!::Function, evaluator!::Function; temperature::Float64 = 1.0, n_steps::Int64 = 0) = MonteCarloDriver(sampler!, evaluator!, temperature, n_steps)

#TODO: Document function
function run!(state::Common.State, driver::MonteCarloDriver, callbacks::Common.CallbackObject...)
    
    step = 0
    xyz0 = copy(state.xyz)
    ene0 = driver.evaluator!(state, false)
    acceptance_count = 0

    while step < driver.n_steps
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

        @Common.cbcall callbacks step state driver (acceptance_count/step)
    end
end

end