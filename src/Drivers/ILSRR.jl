using ProtoSyn.Calculators: EnergyFunction
using ProtoSyn.Mutators: AbstractMutator

Base.@kwdef mutable struct ILSRRState{T <: AbstractFloat} <: DriverState
    step::Int        = 0
    converged::Bool  = false
    completed::Bool  = false
    stalled::Bool    = false
    acceptance_count = 0
    temperature::T   = T(0.0)
end


mutable struct ILSRR <: Driver
    eval!::Union{Function, EnergyFunction}
    jump!::Union{Function, AbstractMutator, Driver}
    inner_driver!::Driver
    callback::Opt{Callback}
    max_steps::Int
    temperature::Function # Takes step, returns temperature at step
end


function (driver::ILSRR)(pose::Pose)
    
    T = eltype(pose.state)
    driver_state = ILSRRState{T}()
    driver_state.temperature = driver.temperature(0)
    
    previous_state  = copy(pose)
    previous_energy = driver.eval!(pose, update_forces = false)
    driver.callback !== nothing && driver.callback(pose, driver_state)
    
    while driver_state.step < driver.max_steps
            
        driver.inner_driver!(pose)
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
            e = pose.state.e[:Total]
            ProtoSyn.recoverfrom!(pose, previous_state)
        end

        driver_state.step += 1
        driver.callback !== nothing && driver.callback(pose, driver_state)

        # Jump
        driver.jump!(pose)
    end

    driver_state.completed = true
    return pose
end

function Base.show(io::IO, drv::ILSRR)
    println(" ILSRR Driver")
    println("\nEnergy function : $(drv.eval!)")
    println("\nJump :")
    println("$(drv.jump!)")
    println("\n+----------------------------------------------------------+")
    println("\nInner driver :")
    println("$(drv.inner_driver!)")
    println("\n+----------------------------------------------------------+\n")
    println("$(drv.callback)")
    println(" Temperature:")
    println("$(drv.temperature)")
    println("\n Settings:")
    println("  Max steps: $(drv.max_steps)")
end