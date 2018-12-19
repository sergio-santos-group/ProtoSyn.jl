module MonteCarlo

using ..Aux
using ..Common
using ..Drivers
using Printf

@doc raw"""
    Driver(sampler!::Function, evaluator!::Function, [, temperature::Float64 = 1.0, n_steps::Int64 = 0, callbacks::Tuple{Common.CallbackObject}...])

Define the runtime parameters for the Monte Carlo simulation.
No `sampler!` movement is performed by default, since n_steps = 0.

# Arguments
- `sampler!::Function`: Responsible for generating a new structure to be evaluated. This function should have the following signature:
```
sampler!(state::Common.State)
```
- `evaluator!::Function`: Responsible for evaluating the system energy. This function should have the following signature:
```
evaluator!(state::Common.State, do_forces::Bool)
```
- `temperature::Float64`: (Optional) Temperature of the system, determines acceptance in the Metropolis algorithm (Default: 1.0).
- `n_steps`: (Optional) Total amount of steps to be performed (Default: 0).
- `callbacks`: (Optional) Tuple of [`CallbackObject`](@ref Common)s.

# Examples
```julia-repl
julia> Drivers.MonteCarlo.Driver(my_sampler!, my_evaluator!, 10.0, 1000)
MonteCarlo.Driver(sampler=my_sampler!, evaluator=my_evaluator!, temperature=10.0, n_steps=1000)

julia> Drivers.MonteCarlo.Driver(my_sampler!, my_evaluator!)
MonteCarlo.Driver(sampler=my_sampler!, evaluator=my_evaluator!, temperature=1.0, n_steps=0)
```
!!! tip
    Both `my_sampler!` and `my_evaluator!` functions often contain pre-defined function avaliable in [Mutators](@ref Mutators) and [Forcefield](@ref Forcefield) modules, respectively.

See also: [`run!`](@ref)
"""
mutable struct Driver <: Drivers.AbstractDriver

    run!::Function
    sampler!::Function
    evaluator!::Function
    temperature::Float64
    n_steps::Int64
    evaluate_slope_every::Int64
    evaluate_slope_threshold::Float64
    verbose::Bool
    callbacks::Tuple

end
function Driver(sampler!::Function, evaluator!::Function, temperature::Float64, n_steps::Int64, evaluate_slope_every::Int64, evaluate_slope_threshold::Float64, verbose::Bool, callbacks::Common.CallbackObject...) 
    return Driver(run!, sampler!, evaluator!, temperature, n_steps, evaluate_slope_every, evaluate_slope_threshold, verbose, callbacks)
end
Base.show(io::IO, b::Driver) = print(io, "MonteCarlo.Driver(sampler=$(string(b.sampler!)) evaluator=$(string(b.evaluator!)), temperature=$(b.temperature), n_steps=$(b.n_steps), evaluate_slope_every=$(b.evaluate_slope_every), evaluate_slope_threshold=$(b.evaluate_slope_threshold), verbose=$(b.verbose))")


@doc raw"""
    run!(state::Common.State, driver::Driver[, callbacks::Tuple{Common.CallbackObject}...])

Run the main body of the driver. Creates a new conformation based on `driver.sampler!`, evaluates the new conformation energy using `driver.evaluator!`,
accepting it or not depending on the `driver.temperature` in a Metropolis algorithm. This Monte Carlo process is repeated for `driver.n_steps`, saving the
accepted structures to `state` and calling all the `callbacks`. 

# Arguments
- `state::Common.State`: Current state of the system to be modified.
- `driver::Driver`: Defines the parameters for the MonteCarlo simulation. See [`Driver`](@ref).
- `callbacks::Vararg{Common.CallbackObject, N}`: (Optional) Tuple of [`CallbackObject`](@ref Common)s. If any callbacks pre-exist in the driver, these are added.

The [`CallbackObject`](@ref Common) in this Driver returns the following extra Varargs (in order):
- `acceptance_ratio::Float64`: The acceptance ratio of the simulation, so far, calculated as `number_of_accepted_steps / current_step`.

# Examples
```julia-repl
julia> Drivers.MonteCarlo.run!(state, driver, my_callback1, my_callback2, my_callback3)
```
"""
function run!(state::Common.State, driver::Driver, callbacks::Common.CallbackObject...)
    
    step = 1
    backup = deepcopy(state)
    driver.evaluator!(state, false)
    acceptance_count = 0
    history_x = Vector{Int64}()
    history_y = Vector{Float64}()

    @Common.cbcall driver.callbacks..., callbacks... 0 state driver (acceptance_count/step)
    while step <= driver.n_steps
        driver.sampler!(state)
        driver.evaluator!(state, false)
        
        if (state.energy.eTotal < backup.energy.eTotal) || (rand() < exp(-(state.energy.eTotal - backup.energy.eTotal) / driver.temperature))
            backup = deepcopy(state)
            push!(history_x, step)
            push!(history_y, state.energy.eTotal)
            acceptance_count += 1
        else
            state = deepcopy(backup)
        end
        
        @Common.cbcall driver.callbacks..., callbacks... step state driver (acceptance_count/step)
        
        # Evaluate slope
        if driver.evaluate_slope_every > 1 && length(history_x) > 0 && length(history_x) % driver.evaluate_slope_every == 0
            b::Float64 = Aux.linreg(history_x, history_y)
            if b >= driver.evaluate_slope_threshold
                if driver.verbose
                    printstyled(@sprintf("(%5s) %12d | Slope analysis: %6.3f ▶️ Exiting inner search ✖\n", "MC", step, b), color = :red)
                end
                break
            end
            history_x = Vector{Int64}()
            history_y = Vector{Float64}()
            if driver.verbose
                printstyled(@sprintf("(%5s) %12d | Slope analysis: %6.3f ▶️ Continuing inner search ✔\n", "MC", step, b), color = :green)
            end
        end
        step += 1
    end
end

end