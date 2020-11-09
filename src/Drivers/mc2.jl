using ProtoSyn.Calculators: EnergyFunction
using ProtoSyn.Mutators: AbstractMutator

Base.@kwdef mutable struct MonteCarloState <: DriverState
    step::Int        = 0
    converged::Bool  = false
    completed::Bool  = false
    stalled::Bool    = false
    acceptance_count = 0
end


Base.@kwdef mutable struct MonteCarlo{T<:AbstractFloat} <: Driver
    eval!::Union{Function, EnergyFunction}
    sample!::Union{Function, AbstractMutator}
    max_steps::Int = 0
    temperature::T
end


function (driver::MonteCarlo{T})(cb::Opt{F}, pose::Pose) where {T, F<:Function}
    
    driver_state = MonteCarloState()
    
    previous_state  = copy(pose)
    previous_energy = driver.eval!(pose, update_forces = false)
    cb !== nothing && cb(pose, driver_state)
    
    while driver_state.step < driver.max_steps
            
        driver.sample!(pose)
        sync!(pose)
        energy = driver.eval!(pose, update_forces = false)
        
        n = rand()
        m = exp((-(energy - previous_energy)) / driver.temperature)
        # println("E: $energy | O: $previous_energy | Difference: $(energy - previous_energy) | $n < $m ? $(n < m)")
        if (energy < previous_energy) || (n < m)
            previous_energy = energy
            previous_state = copy(pose)
            driver_state.acceptance_count += 1
        else
            pose = copy(previous_state)
        end

        driver_state.step += 1
        cb !== nothing && cb(pose, driver_state) 
    end
    
    driver_state.completed = true
    driver_state
end

(driver::MonteCarlo)(pose::Pose) = driver(nothing, pose)