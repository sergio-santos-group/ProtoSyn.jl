module MonteCarlo

using ..Aux
using ..Common
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
Base.@kwdef mutable struct DriverConfig{F <: Function, G <: Function, H <: Function} <: Drivers.AbstractDriver
    sampler!::F
    evaluator!::G
    anneal_fcn::H
    n_steps::Int = 0
end

function DriverConfig(sampler!::F, evaluator!::G, temperature::Float64, n_steps::Int64) where {F <: Function, G <: Function}
    DriverConfig(sampler!, evaluator!, (n::Int64)->temperature, n_steps)
end
# Base.show(io::IO, b::DriverConfig) = print(io, "MonteCarlo.DriverConfig(sampler=$(string(b.sampler!)), evaluator=$(string(b.evaluator!)), n_steps=$(b.n_steps), anneal_fcn=$(string(b.anneal_fcn))")
function Base.show(io::IO, b::DriverConfig)
    print(io, "MonteCarlo.DriverConfig")
    for p in fieldnames(DriverState)
        print(io, "\n   $(String(p)) = $(getproperty(b,p))")
    end
end

# TO DO: Documentation
Base.@kwdef mutable struct DriverState <: Drivers.AbstractDriverState
    step::Int64          = 0
    ac_count::Int        = -1.0
    temperature::Float64 = -1.0
    completed::Bool      = false
end

function Base.show(io::IO, b::DriverState)
    print(io, "MonteCarlo.DriverState")
    for p in fieldnames(DriverState)
        print(io, "\n   $(String(p)) = $(getproperty(b,p))")
    end
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
function run!(state::Common.State, driver_config::DriverConfig, callbacks::Common.CallbackObject...)
    
    # Evaluate initial energy and forces
    #   start by calculating nonbonded lists:
    #   a negative cutoff implies all pairwise
    #   interactions are requested.
    if state.nblist != nothing
        state.nblist.cutoff = -1.0
    end
    update_nblist(state,top)
    energy = driver_config.evaluator!(state, false)
    
    # instantiate a new DriverState object.
    # By default, no optimization step has yet been taken
    # apart from calculating the energy and forces for the
    # input state
    driver_state = DriverState()
    
    # create a copy of the input state
    prev_state = Common.State(state.size)
    
    # initialize auxilliary variables
    driver_state.ac_count = 0    # accepted counter
    
    # call "callback" functions
    @Common.cbcall callbacks state driver_state driver_config

    #region MAINLOOP
    while driver_state.step < driver_config.n_steps
        driver_state.step += 1
        
        # sample new configuration
        driver_config.sampler!(state)

        # evaluate energy of new configuration
        energy = driver_config.evaluator!(state, false)

        # calculate temperature for current step
        driver_state.temperature = driver_config.anneal_fcn(driver_state.step)

        ΔE = energy - prev_state.energy.total
        if (ΔE < 0.0) || (rand() < exp(-ΔE/driver_state.temperature) )
            # since the new configuration was accepted,
            # copy it to the prev_state and increment
            # the accepted counter
            copy!(prev_state.xyz, state.xyz)
            copy!(prev_state.energy, state.energy)
            # copy!(prev_state, state)
            driver_state.ac_count += 1
        else
            # otherwise, revert to the previous state
            copy!(state.xyz, prev_state.xyz)
            copy!(state.energy, prev_state.energy)
            # copy!(state, prev_state)
        end
        
        # update driver state and call calback functions (if any)
        @Common.cbcall callbacks state driver_state driver_config

    end
    #endregion

    # update driver state and return
    driver_state.completed = true

    return  driver_state
    
end

end