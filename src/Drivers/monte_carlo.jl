using ProtoSyn.Calculators: EnergyFunction
using ProtoSyn.Mutators: AbstractMutator

Base.@kwdef mutable struct MonteCarloState{T <: AbstractFloat} <: DriverState
    step::Int        = 0
    converged::Bool  = false
    completed::Bool  = false
    stalled::Bool    = false
    acceptance_count = 0
    temperature::T   = 0.0
end


mutable struct MonteCarlo <: Driver
    eval!::Union{Function, EnergyFunction}
    sample!::Union{Function, AbstractMutator, Driver}
    callback::Opt{Callback}
    max_steps::Int
    temperature::Function
end


function (driver::MonteCarlo)(pose::Pose)
    
    T = eltype(pose.state)
    driver_state = MonteCarloState{T}()
    driver_state.temperature = driver.temperature(0)
    
    previous_state  = copy(pose)
    previous_energy = driver.eval!(pose, update_forces = false)
    driver.callback !== nothing && driver.callback(pose, driver_state)
    
    while driver_state.step < driver.max_steps
            
        driver.sample!(pose)
        sync!(pose)
        energy = driver.eval!(pose, update_forces = false)
        
        n = rand()
        driver_state.temperature = driver.temperature(driver_state.step)
        m = exp((-(energy - previous_energy)) / driver_state.temperature)
        if (energy < previous_energy) || (n < m)
            previous_energy = energy
            previous_state = copy(pose)
            driver_state.acceptance_count += 1
        else
            pose = copy(previous_state)
        end

        driver_state.step += 1
        driver.callback !== nothing && driver.callback(pose, driver_state) 
    end
    
    driver_state.completed = true
    driver_state
end