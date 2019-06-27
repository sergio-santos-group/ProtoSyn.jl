module ILSRR

using ..Aux
using ..Common
using ..Drivers
using ..Abstract

#= ----------------------------------------------------------------------------
ILSRR (Iterated Local Search with Random Restarts)

# Algorithm explanation:
The ILSRR algorithm searchs the local conformational neighbourhood using
an "inner driver", such as Monte Carlo or even Steepest Descent. The
objective of the inner driver is to find local optimums. The encompassing
outer cycle is responsible to perform "perturbations": big mutations in the
structure or nature of the molecule being sampled in order to "kick" a solution
out of the local optimum, essentially overcomming energy barriers.
The starting point of the inner cycle is called the "home base", and is re-defined
whenever a new optimum is found or the energy gap in comparison to the last
home base is small enough to pass the Metropolis criterium.

# Convergence criteria
1. Number of `n_steps` performed;
2. Number of `n_stalls` cycles failed to produce a new Home Base;
---------------------------------------------------------------------------- =#


@doc raw"""
    DriverConfig(inner_cycle_config::Abstract.DriverConfig, pertubator::Abstract.Sampler[, temperature::Union{Function, Float64} = 0.0, n_steps::Int64 = 0, stall_limit::Int64 = 0, callbacks::Vector{Abstract.CallbackObject} = []])

Define the runtime parameters for the ILSRR algorithm.

# Arguments
- `inner_cycle_config::Abstract.DriverConfig`: Responsible for driving the inner cycle of the ILSRR algorithm. Should be a Driver, such as [`MonteCarlo`](@ref Drivers)
- `pertubator::Abstract.Sampler`: Responsible for performing "big" conformational changes in the system.
- `temperature::Union{Function, Float64}`: (Optional) Temperature for the Metropolis criteria.
This parameter is always a function who should follow the following signature:
```
function (n::Int64)::Float64
```
When a Float64 is passed, by default, ProtoSyn will construct an unnamed function who always returns that Float64 number every step (Default: 0.0)
- `n_steps`: (Optional) Total amount of outer cycles to be performed (Default: 0).
- `stall_limit`: (Optional) If defined, will return after `stall_limit` jumps who consecutively fail to produce a new optimum (Default: 0).
- `callbacks`: (Optional) Vector of [`CallbackObject`](@ref Common)s.

# Examples
```julia-repl
julia> Drivers.ILSRR.Driver(inner_cycle_driver, my_pertubator)
```
See also: [`run!`](@ref)
"""
mutable struct DriverConfig{F <: Function, T <: Abstract.CallbackObject} <: Abstract.DriverConfig

    inner_driver_config::Abstract.DriverConfig # Required
    perturbator::Abstract.Sampler              # Required
    anneal_fcn::F                              # Default: constant temperature 0.0
    n_steps::Int64                             # Default = 0
    stall_limit::Int64                         # Default = 0
    callbacks::Vector{T}

    DriverConfig(; inner_driver_config::Abstract.DriverConfig,
        perturbator::Abstract.Sampler,
        temperature::Union{Float64, Function} = 0.0,
        n_steps::Int64 = 0,
        stall_limit::Int64 = 0,
        callbacks::Vector{<:Abstract.CallbackObject} = Vector{Common.CallbackObject}()) = begin

        if typeof(temperature) == Float64
            new{Function, Abstract.CallbackObject}(inner_driver_config,
                perturbator,
                function constant_temperature(n::Int64) temperature end,
                n_steps,
                stall_limit,
                callbacks)
        else
            new{Function, Abstract.CallbackObject}(inner_driver_config,
                perturbator,
                temperature,
                n_steps,
                stall_limit,
                callbacks)
        end
    end
end


@doc raw"""
    DriverState(; best_state::Union{Common.State, Nothing} = nothing, home_state::Union{Common.State, Nothing} = nothing, step::Int64 = 0, n_stalls::Int64 = 0, temperature::Float64 = -1.0, completed::Bool = false, stalled::Bool = false)

Define the runtime parameters given by the Monte Carlo simulation.

# Arguments
- `best_state::Union{Common.State, Nothing}`: Current best state of simulation.
- `home_state::Union{Common.State, Nothing}`: Current home state of simulation.
- `step::Int64`: Current step of the simulation.
- `n_stalls::Int64`: Current number of consecutive perturbations who failed to produce a new optimum. 
- `temperature::Float64`: Current temperature of the system.
- `completed::Bool`: Current completness status of the system.
- `stalled::Bool`: Current stalled status of the system.

# Examples
```julia-repl
julia> Drivers.ILSRR.DriverState()
ProtoSyn.Drivers.ILSRR.DriverState
   best_state = nothing
   home_state = nothing
   step = 0
   n_stalls = 0
   temperature = -1.0
   completed = false
   stalled = false
```
See also: [`run!`](@ref)
"""
Base.@kwdef mutable struct DriverState <: Abstract.DriverState

    # Parameter:                             # Default value:
    best_state::Union{Common.State, Nothing} = nothing
    home_state::Union{Common.State, Nothing} = nothing
    step::Int                                = 0
    n_stalls::Int                            = 0
    temperature::Float64                     = -1.0
    completed::Bool                          = false
    stalled::Bool                            = false
end


# ----------------------------------------------------------------------------------------------------------
#                                                   RUN

@doc raw"""
    run!(state::Common.State, driver_config::DriverConfig)

Run the main body of the Driver.

# Arguments
- `state::Common.State`: Current state of the system to be modified.
- `driver_config::DriverConfig`: Defines the parameters for the ILSRR algorithm. See [`DriverConfig`](@ref).

# Examples
```julia-repl
julia> Drivers.ILSRR.run!(state, ilsrr_driver)
```
"""
function run!(state::Common.State, driver_config::DriverConfig)

    driver_state = DriverState()
    
    inner_driver_config = driver_config.inner_driver_config
    
    let n_steps=inner_driver_config.n_steps
        inner_driver_config.n_steps = 0
        Drivers.run!(state, inner_driver_config)
        inner_driver_config.n_steps = n_steps
    end
    
    driver_state.best_state = Common.State(state)
    driver_state.home_state = Common.State(state)
    driver_state.completed = driver_state.step == driver_config.n_steps

    Common.@cbcall driver_config.callbacks state driver_state
    
    R = 0.0083144598 # kJ mol-1 K-1

    #region MAINLOOP
    while !(driver_state.completed || driver_state.stalled)
        
        # this driver should make multiple small tweaks
        # to the state
        typeof(inner_driver_config).name.module.run!(state, inner_driver_config)
        
        driver_state.step += 1
        driver_state.temperature = driver_config.anneal_fcn(driver_state.step)
        
        if state.energy.total < driver_state.best_state.energy.total
            # save this state as the best and make
            # it the new homebase
            Common.@copy driver_state.best_state state energy xyz
            Common.@copy driver_state.home_state state energy xyz
            driver_state.n_stalls = 0
        else
            # otherwise, a new homebase may be created according
            # to the Metropolis criterium
            β = driver_state.temperature != 0.0 ? 1/(R * driver_state.temperature) : Inf
            ΔE = state.energy.total - driver_state.home_state.energy.total
            if (ΔE <= 0.0) || (rand() < exp(-ΔE*β) )
                Common.@copy driver_state.home_state state energy xyz
                driver_state.n_stalls = 0
            else
                # if the criterium was not accepted, revert
                # to the homebase
                Common.@copy state driver_state.home_state energy xyz
                driver_state.n_stalls += 1
            end
        end
        
        driver_state.stalled = driver_state.n_stalls == driver_config.stall_limit
        driver_state.completed = driver_state.step == driver_config.n_steps
        
        # make a large perturbation to the state
        if !(driver_state.completed || driver_state.stalled)
            driver_config.perturbator.apply!(state, driver_config.perturbator.mutators)
        end

        Common.@cbcall driver_config.callbacks state driver_state
    end
    #endregion

    # before returning, save the best state
    if driver_state.best_state.energy.total < state.energy.total
        Common.@copy state driver_state.best_state energy xyz
    end

    return driver_state
end # end function

end # end module