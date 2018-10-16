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

#TODO: Document function
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