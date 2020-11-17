using ProtoSyn.Calculators: EnergyFunction
using ProtoSyn.Mutators: AbstractMutator

Base.@kwdef mutable struct MonteCarloState <: DriverState
    step::Int        = 0
    converged::Bool  = false
    completed::Bool  = false
    stalled::Bool    = false
    acceptance_count = 0
end


mutable struct MonteCarlo{T<:AbstractFloat} <: Driver
    eval!::Union{Function, EnergyFunction}
    sample!::Union{Function, AbstractMutator, Driver}
    callback::Opt{Callback}
    max_steps::Int
    temperature::T
end


function (driver::MonteCarlo{T})(pose::Pose) where {T}
    
    driver_state = MonteCarloState()
    
    previous_state  = copy(pose)
    previous_energy = driver.eval!(pose, update_forces = false)
    driver.callback !== nothing && driver.callback.event(pose, driver_state)
    
    while driver_state.step < driver.max_steps
            
        driver.sample!(pose)
        sync!(pose)
        energy = driver.eval!(pose, update_forces = false)
        
        n = rand()
        m = exp((-(energy - previous_energy)) / driver.temperature)
        if (energy < previous_energy) || (n < m)
            previous_energy = energy
            previous_state = copy(pose)
            driver_state.acceptance_count += 1
        else
            pose = copy(previous_state)
        end

        driver_state.step += 1
        driver.callback !== nothing && driver.callback.event(pose, driver_state) 
    end
    
    driver_state.completed = true
    driver_state
end