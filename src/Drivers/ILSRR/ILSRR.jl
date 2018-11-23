module ILSRR

using ..Aux
using ..Common
using ..Drivers
using Printf

mutable struct Driver <: Drivers.AbstractDriver

    inner_cycle_driver::Drivers.AbstractDriver
    evaluator!::Function
    perturbator!::Function
    temperature::Float64
    n_steps::Int64
    run!::Function

end
Driver(inner_cycle_driver::Drivers.AbstractDriver, evaluator!::Function, perturbator!::Function; n_outer_steps::Int64 = 2, outer_temperature::Float64 = 0.0) = Driver(inner_cycle_driver, evaluator!, perturbator!, n_outer_steps, outer_temperature, run!)



function run!(state::Common.State, driver::Driver, callbacks::Common.CallbackObject...)

    save_inner_best = @Common.callback 1 function cb_save(step::Int64, st::Common.State, dr::Drivers.AbstractDriver, args...)
        if st.energy.eTotal < inner_best.energy.eTotal
            inner_best = deepcopy(st)
        end
    end

    step = 0
    driver.evaluator!(state, false)
    inner_best = deepcopy(state)
    outer_best = deepcopy(state)
    homebase   = deepcopy(state)
    while step < driver.n_steps
        driver.inner_cycle_driver.run!(state, driver.inner_cycle_driver, save_inner_best, callbacks...)
        state = deepcopy(inner_best)

        if state.energy.etotal < outer_best.energy.eTotal
            outer_best = deepcopy(state)
        end

        if state.energy.eTotal < homebase.energy.eTotal || (rand() < exp(-(state.energy.eTotal - homebase.energy.eTotal) / driver.temperature))
            homebase = deepcopy(state)
        else
            state = deepcopy(homebase)
        end

        if step != n_steps
            movs = driver.perturbator!(state)
            driver.evaluator!(state, false)
            inner_best = deepcopy(state)
        end

        @Common.cbcall callbacks step state driver movs
        outer_step += 1
    end
end

end