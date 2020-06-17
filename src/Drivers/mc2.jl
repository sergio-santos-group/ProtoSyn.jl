

#------------------------------------------------
Base.@kwdef mutable struct MonteCarloState <: IDriverState
    step::Int       = 0
    converged::Bool = false
    completed::Bool = false
    stalled::Bool   = false

    #stepsize::Float64 =  1.0
    #max_force::Float64 = -1.0
    acceptance_count = 0
end

#------------------------------------------------
Base.@kwdef mutable struct MonteCarlo{T<:AbstractFloat} <: IDriver
    eval!::Function
    sample!::Function
    max_steps::Int = 0
    #pairlist_freq::Int = 0

    #force_tolerance::Float64 = 1.0
    #max_displacement::Float64 = 0.1
    temperature::T
end


# function (driver::MonteCarlo{T})(cb::Opt{F}, state::Calculators.State{T}) where {T, F<:Function}
# function (driver::MonteCarlo{T})(cb::Opt{F}, state::State{T}) where {T, F<:Function}
    
#     driver_state = MonteCarloState()
    
#     coords0 = copy(state.x)
#     energy0 = driver.eval!(state, false)

#     cb!==nothing && cb(state, driver_state)
    
#     while driver_state.step < driver.max_steps
            
#         driver.sample!(state)
#         energy1 = driver.eval!(state, false)
        
#         if (energy1 < energy0) || (rand() < exp((energy0-energy1)/driver.temperature))
#             energy0 = energy1
#             copyto!(coords0, state.x)
#             driver_state.acceptance_count += 1
#         else
#             copyto!(state.x, coords0)
#         end

#         driver_state.step += 1
#         cb!==nothing && cb(state, driver_state)
    
#     end
    
#     driver_state.completed = true
#     driver_state
# end

(driver::MonteCarlo)(s::State) = driver(nothing, s)

function (driver::MonteCarlo{T})(cb::Opt{F}, state::AbstractRepresentation) where {T, F<:Function}
    
    driver_state = MonteCarloState()
    
    prepare(state)
    prevstate = save(state)
    prevenergy = driver.eval!(state, false)

    cb!==nothing && cb(state, driver_state)
    
    while driver_state.step < driver.max_steps
            
        driver.sample!(state)
        energy = driver.eval!(state, false)
        
        if (energy < prevenergy) || (rand() < exp((prevenergy-energy)/driver.temperature))
            energy0 = energy
            save!(prevstate, state)
            driver_state.acceptance_count += 1
        else
            restore!(state, prevstate)
        end

        driver_state.step += 1
        cb!==nothing && cb(state, driver_state)
    
    end
    
    driver_state.completed = true
    driver_state
end
