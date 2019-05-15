module ILSRR
const id = :ILSRR

using ..Aux
using ..Common
using ..Drivers
using ..Abstract

#= ----------------------------------------------------------------------------
ILSRR Example
(Iterated Local Search with Random Restarts)

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
1. Number of `n_steps` performed
2. Number of `n_stalls` cycles failed to produce a new Home Base
---------------------------------------------------------------------------- =#


@doc raw"""
    Driver(inner_cycle_driver!::Drivers.AbstractDriver, evaluator!::Function, pertubator!::Function[, temperature::Float64 = 0.0, n_steps::Int64 = 0, callbacks::Tuple{Common.CallbackObject}...])

Define the runtime parameters for the ILSRR algorithm.

# Arguments
- `inner_cycle_driver!::Driver.AbstractDriver`: Responsible for driving the inner cycle of the ILSRR algorithm. Should be a Driver, such as [`MonteCarlo`](@ref Drivers)
- `evaluator!::Function`: Responsible for evaluating the current `state.energy`. This function should have the following signature:
```
evaluator!(state::Common.State, do_forces::Bool)
```
- `pertubator!::Function`: Responsible for performing conformational changes in the system. It's usually an aggregation of [Mutators](@ref Mutators).
- `temperature::Float64`: (Optional) Temperature for the Metropolis criteria when performing system perturbation (Default: 0.0).
- `n_steps`: (Optional) Total amount of outer cycles to be performed (Default: 0).
- `continue_after_n_attemps`: (Optional) If defined, will reset to initial structure after `continue_after_n_attemps` jumps who consecutively fail to produce a new optimum (Default: 0).
- `callbacks`: (Optional) Tuple of [`CallbackObject`](@ref Common)s.

# Examples
```julia-repl
julia> Drivers.ILSRR.Driver(inner_cycle_driver, my_evaluator!, my_pertubator)
ILSRR.Driver(evaluator=my_evaluator!, n_steps=100, f_tol=1e-3, max_step=0.1)

julia> Drivers.ILSRR.Driver(inner_cycle_driver, my_evaluator!, my_pertubator, 300.0, 10)
ILSRR.Driver(evaluator=my_evaluator!, temperature=300.0, n_steps=10)
```
!!! tip
    The `my_evaluator!` function often contains an aggregation of pre-defined functions avaliable in [Forcefield](@ref Forcefield). It is possible to combine such functions using the [`@faggregator`](@ref Common) macro.

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

#TO DO: Documentation
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
    run!(state::Common.State, driver::SteepestDescentDriver[, callback::Union{Common.CallbackObject, Nothing} = nothing])

Run the main body of the Driver.

# Arguments
- `state::Common.State`: Current state of the system to be modified.
- `driver::SteepestDescentDriver`: Defines the parameters for the ILSRR algorithm. See [`Driver`](@ref).
- `callbacks::Vararg{Common.CallbackObject, N}`: (Optional) Tuple of [`CallbackObject`](@ref Common)s.

!!! tip
    The callback function often contains a [Print](@ref) function.

# Examples
```julia-repl
julia> Drivers.ILSRR.run(state, ilsrr_driver, callback1, callback2, callback3)
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
            n_stalls = 0
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