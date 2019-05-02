module MonteCarlo

using ..Aux
using ..Common
using ..Drivers
using Printf

@doc raw"""
    Driver(sampler!::Function, evaluator!::Function, [, temperature::Float64 = 1.0, n_steps::Int64 = 0, callbacks::Tuple{Common.CallbackObject}...])

Define the runtime parameters for the Monte Carlo simulation.
No `sampler!` movement is performed by default, since n_steps = 0.

# Arguments
- `sampler!::Function`: Responsible for generating a new structure to be evaluated. This function should have the following signature:
```
sampler!(state::Common.State)
```
- `evaluator!::Function`: Responsible for evaluating the system energy. This function should have the following signature:
```
evaluator!(state::Common.State, do_forces::Bool)
```
- `temperature::Float64`: (Optional) Temperature of the system, determines acceptance in the Metropolis algorithm (Default: 1.0).
- `n_steps`: (Optional) Total amount of steps to be performed (Default: 0).
- `callbacks`: (Optional) Tuple of [`CallbackObject`](@ref Common)s.

# Examples
```julia-repl
julia> Drivers.MonteCarlo.Driver(my_sampler!, my_evaluator!, 10.0, 1000)
MonteCarlo.Driver(sampler=my_sampler!, evaluator=my_evaluator!, temperature=10.0, n_steps=1000)

julia> Drivers.MonteCarlo.Driver(my_sampler!, my_evaluator!)
MonteCarlo.Driver(sampler=my_sampler!, evaluator=my_evaluator!, temperature=1.0, n_steps=0)
```
!!! tip
    Both `my_sampler!` and `my_evaluator!` functions often contain pre-defined function avaliable in [Mutators](@ref Mutators) and [Forcefield](@ref Forcefield) modules, respectively.

See also: [`run!`](@ref)
"""
Base.@kwdef mutable struct DriverConfig{F <: Function, G <: Function, H <: Function} <: Drivers.AbstractDriver
    sampler!::F
    evaluator!::G
    anneal_fcn::H
    n_steps::Int = 0
end

function DriverConfig(sampler!::F, evaluator!::G, temperature::Float64, n_steps::Int64) where {F <: Function, G <: Function}
    DriverConfig(sampler!, evaluator!, (n::Int64)->temperature, n_steps)
end
Base.show(io::IO, b::DriverConfig) = print(io, "MonteCarlo.DriverConfig(sampler=$(string(b.sampler!)), evaluator=$(string(b.evaluator!)), anneal_fcn=$(string(b.anneal_fcn)), n_steps=$(b.n_steps)")

# TO DO: Documentation
Base.@kwdef mutable struct DriverState
    step::Int64          = 0
    ac_ratio::Float64    = 0.0
    temperature::Float64 = 0.0
end
Base.show(io::IO, b::DriverState) = print(io, "MonteCarlo.DriverState(step=$(b.step), ac_ratio=$(b.ac_ratio), temperature=$(b.temperature))")

# ----------------------------------------------------------------------------------------------------------
#                                                   RUN

@doc raw"""
    run!(state::Common.State, driver::Driver[, callbacks::Tuple{Common.CallbackObject}...])

Run the main body of the driver. Create a new conformation based on `driver.sampler!`, evaluates the new conformation energy using `driver.evaluator!`,
accepting it or not depending on the `driver.temperature` in a Metropolis algorithm. This Monte Carlo process is repeated for `driver.n_steps`, saving the
accepted structures to `state` and calling all the `callbacks`. 

# Arguments
- `state::Common.State`: Current state of the system to be modified.
- `driver::Driver`: Defines the parameters for the MonteCarlo simulation. See [`Driver`](@ref).
- `callbacks::Vararg{Common.CallbackObject, N}`: (Optional) Tuple of [`CallbackObject`](@ref Common)s. If any callbacks pre-exist in the driver, these are added.

The [`CallbackObject`](@ref Common) in this Driver returns the following extra Varargs (in order):
- `acceptance_ratio::Float64`: The acceptance ratio of the simulation, so far, calculated as `number_of_accepted_steps / current_step`.

# Examples
```julia-repl
julia> Drivers.MonteCarlo.run!(state, driver, my_callback1, my_callback2, my_callback3)
```
"""
function run!(state::Common.State, driver_config::DriverConfig, callbacks::Common.CallbackObject...)
    
    driver_config.evaluator!(state, false)
    
    ac_count = 0
    driver_state = DriverState(1, 0.0, driver_config.anneal_fcn(1))
    backup_state = Common.State(state.size)
    copy!(backup_state, state)

    @Common.cbcall callbacks state driver_config driver_state
    while driver_state.step <= driver_config.n_steps
        driver_config.sampler!(state)
        driver_config.evaluator!(state, false)

        driver_state.temperature = driver_config.anneal_fcn(driver_state.step)

        #@metropolis state.energy.eTotal backup.energy.eTotal temperature
        if (state.energy.eTotal < backup_state.energy.eTotal) || (rand() < exp(-(state.energy.eTotal - backup_state.energy.eTotal) / driver_state.temperature)) # Metropolis
            copy!(backup_state, state)
            ac_count += 1
        else
            copy!(state, backup_state)
        end
        
        driver_state.ac_ratio = ac_count / driver_state.step
        @Common.cbcall callbacks state driver_config driver_state
        driver_state.step += 1
    end
end

end