module ILSRR

using ..Aux
using ..Common
using ..Drivers
using Printf


@doc raw"""
    Driver(inner_cycle_driver!::Drivers.AbstractDriver, evaluator!::Function, pertubator!::Function[, temperature::Float64 = 0.0, n_steps::Int64 = 0, callbacks::Tuple{Common.CallbackObject}...])

Define the runtime parameters for the ILSRR algorithm.

# Arguments
- `inner_cycle_driver!::Driver.AbstractDriver`: Responsible for driving the inner cycle of the ILSRR algorithm. Should be a Driver, such as [`MonteCarlo`](@ref Drivers)
- `evaluator!::Function`: Responsible for evaluating the current `state.energy`. This function should have the following signature:
```
evaluator!(state::Common.State, do_forces::Bool)
```
- `pertubator!::Function`: Responsible for performing conformational changes in the system. It's usually an aggregation of [Mutators](@ref Mutators).
- `temperature::Float64`: (Optional) Temperature for the Metropolis criteria when performing system perturbation (Default: 0.0)
- `n_steps`: (Optional) Total amount of outer cycles to be performed (Default: 0).
- `callbacks`: (Optional) Tuple of [`CallbackObject`](@ref Common)s.

# Examples
```julia-repl
julia> Drivers.ILSRR.Driver(inner_cycle_driver, my_evaluator!, my_pertubator)
ILSRR.Driver(evaluator=my_evaluator!, n_steps=100, f_tol=1e-3, max_step=0.1)

julia> Drivers.ILSRR.Driver(inner_cycle_driver, my_evaluator!, my_pertubator, 300.0, 10)
ILSRR.Driver(evaluator=my_evaluator!, temperature=300.0, n_steps=10)
```
!!! tip
    The `my_evaluator!` function often contains an aggregation of pre-defined functions avaliable in [Forcefield](@ref Forcefield). It is possible to combine such functions using the [`@faggregator`](@ref Common) macro.

See also: [`run!`](@ref)
"""
mutable struct Driver <: Drivers.AbstractDriver

    run!::Function
    inner_cycle_driver::Drivers.AbstractDriver
    evaluator!::Function
    perturbator!::Function
    temperature::Float64
    n_steps::Int64
    callbacks::Tuple

end
Driver(inner_cycle_driver::Drivers.AbstractDriver, evaluator!::Function, perturbator!::Function, temperature::Float64 = 0.0, n_steps::Int64 = 0, callbacks::Common.CallbackObject...) = Driver(run!, inner_cycle_driver, evaluator!, perturbator!, temperature, n_steps, callbacks)
Base.show(io::IO, b::Driver) = print(io, "ILSRR.Driver(inner_cycle_driver=$(b.inner_cycle_driver), evaluator=$(string(b.evaluator!)), n_steps=$(b.n_steps), f_tol=$(b.f_tol), max_step=$(b.max_step))")


@doc raw"""
    run!(state::Common.State, driver::SteepestDescentDriver[, callback::Union{Common.CallbackObject, Nothing} = nothing])

Run the main body of the Driver.

# Arguments
- `state::Common.State`: Current state of the system to be modified.
- `driver::SteepestDescentDriver`: Defines the parameters for the ILSRR algorithm. See [`Driver`](@ref).
- `callbacks::Vararg{Common.CallbackObject, N}`: (Optional) Tuple of [`CallbackObject`](@ref Common)s.

!!! tip
    The callback function often contains a [Print](@ref) function.

# Examples
```julia-repl
julia> Drivers.ILSRR.run(state, ilsrr_driver, callback1, callback2, callback3)
```
"""
function run!(state::Common.State, driver::Driver, callbacks::Common.CallbackObject...)

    save_inner_best = @Common.callback 1 function cb_save(step::Int64, st::Common.State, dr::Drivers.AbstractDriver, args...)
        if st.energy.eTotal < inner_best.energy.eTotal
            inner_best = deepcopy(st)
        end
    end

    step = 1
    driver.evaluator!(state, false)
    inner_best = deepcopy(state)
    homebase   = deepcopy(state)
    # @Common.cbcall driver.callbacks..., callbacks... step state driver
    while step <= driver.n_steps
        println(@sprintf("\n(%5s) %12s \n%s", "ILSRR", @sprintf("Step: %4d", step), "-"^150))

        driver.inner_cycle_driver.run!(state, driver.inner_cycle_driver, save_inner_best)
        state = deepcopy(inner_best)

        step += 1
        if step <= driver.n_steps
            @Common.cbcall driver.callbacks..., callbacks... step state driver
        end

        if state.energy.eTotal < homebase.energy.eTotal || (rand() < exp(-(state.energy.eTotal - homebase.energy.eTotal) / driver.temperature))
            println(@sprintf "(ILSRR) New homebase defined: ⚡E: %10.3e (old) ▶️ %10.3e (new)" homebase.energy.eTotal state.energy.eTotal)
            homebase = deepcopy(state)
        else
            println(@sprintf "(ILSRR) Recovering to previous homebase: ⚡E: %10.3e (actual) ▶️ %10.3e (new)" state.energy.eTotal homebase.energy.eTotal)
            state = deepcopy(homebase)
        end

        if (step-1) != driver.n_steps
            driver.perturbator!(state)
            driver.evaluator!(state, false)
            inner_best = deepcopy(state)
        end
    end
end

end