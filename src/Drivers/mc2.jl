

#------------------------------------------------
Base.@kwdef mutable struct MonteCarloState <: IDriverState
    step::Int64        = 0
    converged::Bool    = false
    completed::Bool    = false
    stalled::Bool      = false

    #stepsize::Float64 =  1.0
    #max_force::Float64 = -1.0
    acceptance_count = 0
end

#------------------------------------------------
Base.@kwdef mutable struct MonteCarlo{F <: Function, G<:Function} <: IDriver
    eval!::F
    sample!::G
    max_steps::Int = 0
    #pairlist_freq::Int = 0

    #force_tolerance::Float64 = 1.0
    #max_displacement::Float64 = 0.1
    temperature::AbstractFloat
end


function (driver::MonteCarlo)(cb::Opt{F}, state::State) where {F <: Function}
    
    driver_state = MonteCarloState()
    
    old_coords = copy(state.coords)
    old_energy = driver.eval!(state, false)

    cb!==nothing && cb(state, driver_state)
    
    while driver_state.step < driver.max_steps
            
        driver.sample!(state)
        new_energy = driver.eval!(state, false)
        
        if (new_energy < old_energy) || (rand() < exp((old_energy-new_energy)/driver.temperature))
            old_energy = new_energy
            copyto!(old_coords, state.coords)
            driver_state.acceptance_count += 1
        else
            copyto!(state.coords, old_coords)
        end

        driver_state.step += 1
        cb!==nothing && cb(state, driver_state)
    
    end
    
    driver_state.completed = true
    driver_state
end

(driver::MonteCarlo)(s::State) = driver(nothing, s)