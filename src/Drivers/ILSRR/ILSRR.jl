module ILSRR

using ..Aux
using ..Common
using ..Drivers
using Printf

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
    @Common.cbcall driver.callbacks..., callbacks... step state driver
    while step <= driver.n_steps
        driver.inner_cycle_driver.run!(state, driver.inner_cycle_driver, save_inner_best)
        state = deepcopy(inner_best)

        if state.energy.eTotal < homebase.energy.eTotal || (rand() < exp(-(state.energy.eTotal - homebase.energy.eTotal) / driver.temperature))
            homebase = deepcopy(state)
        else
            state = deepcopy(homebase)
        end

        if step != driver.n_steps
            driver.perturbator!(state)
            driver.evaluator!(state, false)
            inner_best = deepcopy(state)
        end

        step += 1
        if step <= driver.n_steps
            @Common.cbcall driver.callbacks..., callbacks... step state driver
        end
    end
end

end