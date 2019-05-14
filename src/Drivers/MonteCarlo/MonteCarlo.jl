module MonteCarlo
const id = :MC

using ..Aux
using ..Common
using ..Abstract
using ..Drivers

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
mutable struct DriverConfig{F <: Function, T <: Abstract.CallbackObject} <: Abstract.DriverConfig
    sampler::Abstract.Sampler                # Required
    evaluator::Abstract.Evaluator            # Required
    anneal_fcn::F                            # Default: constant temperature 0.0
    n_steps::Int                             # Default: 0
    callbacks::Vector{T} # Default: empty

    DriverConfig(; sampler::Abstract.Sampler,
        evaluator::Abstract.Evaluator,
        temperature::Union{Float64, Function} = 0.0,
        n_steps::Int64 = 0,
        callbacks::Vector{<:Abstract.CallbackObject} = Vector{Common.CallbackObject}()) = begin

        if typeof(temperature) == Float64
            new{Function, Abstract.CallbackObject}(sampler,
                evaluator,
                function constant_temperature(n::Int64) temperature end,
                n_steps,
                callbacks)
        else
            new{Function, Abstract.CallbackObject}(sampler,
                evaluator,
                temperature,
                n_steps,
                callbacks)
        end
    end
end


# TODO: Documentation
Base.@kwdef mutable struct DriverState <: Abstract.DriverState
    
    # Parameter:         # Default:
    step::Int64          = 0
    ac_count::Int        = -1
    temperature::Float64 = -1.0
    completed::Bool      = false
end


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
function run!(state::Common.State, driver_config::DriverConfig)
    
    # Evaluate initial energy and forces
    #   start by calculating nonbonded lists:
    #   a negative cut_off implies all pairwise
    #   interactions are requested.
    if state.nblist != nothing
        state.nblist.cut_off = -1.0
    end
    Common.update_nblist!(state)
    energy = driver_config.evaluator.evaluate!(state, driver_config.evaluator.components, false)
    
    # instantiate a new DriverState object.
    # By default, no optimization step has yet been taken
    # apart from calculating the energy and forces for the
    # input state
    driver_state = DriverState()
    
    # create a copy of the input state
    prev_state = Common.State(state.size)
    copy!(prev_state.xyz, state.xyz)
    copy!(prev_state.energy, state.energy)
    
    # initialize auxilliary variables
    driver_state.ac_count = 0    # accepted counter
    
    # call "callback" functions
    Common.@cbcall driver_config.callbacks state driver_state

    R = 0.0083144598 # kJ mol-1 K-1


    #region MAINLOOP
    while driver_state.step < driver_config.n_steps
        driver_state.step += 1
        
        # sample new configuration
        driver_config.sampler.apply!(state, driver_config.sampler.mutators)
        
        # evaluate energy of new configuration
        energy = driver_config.evaluator.evaluate!(state, driver_config.evaluator.components, false)
        
        # calculate temperature for current step
        # if T ∈ ]0, +∞[ : exp(-ΔE/(R*T)) ∈ ]0, 1]                     -> Metropolis MAY be accepted
        # if T ∈ ]-∞, 0[ : exp(-ΔE/(R*T)) ∈ [1, +∞[                    -> Metropolis is ALWAYS accepted
        # if T == 0 : exp(-ΔE/(R*T)) = 0 because one takes 1/0 as Inf  -> Metropolis is NEVER accepted
        #   despite not being possible division by zero. 
        driver_state.temperature = driver_config.anneal_fcn(driver_state.step)
        β = driver_state.temperature != 0.0 ? 1/(R * driver_state.temperature) : Inf

        ΔE = energy - prev_state.energy.total
        if (ΔE <= 0.0) || (rand() < exp(-ΔE*β) )
            # since the new configuration was accepted,
            # copy it to the prev_state and increment
            # the accepted counter
            copy!(prev_state.xyz, state.xyz)
            copy!(prev_state.energy, state.energy)
            driver_state.ac_count += 1
        else
            # otherwise, revert to the previous state
            copy!(state.xyz, prev_state.xyz)
            copy!(state.energy, prev_state.energy)
        end
        
        # update driver state and call calback functions (if any)
        Common.@cbcall driver_config.callbacks state driver_state

    end
    #endregion

    # update driver state and return
    driver_state.completed = true

    return  driver_state
    
end

end