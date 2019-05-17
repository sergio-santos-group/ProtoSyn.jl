module MonteCarlo
const id = :MC

using ..Aux
using ..Common
using ..Abstract
using ..Drivers

#= ----------------------------------------------------------------------------
Monte Carlo

# Algorithm explanation:
The Monte Carlo algorithm samples the local conformational neighbourhood.
It is usually comprised of two indispensable components: a sampler and an
evaluator. The sampler is responsible to proposing a new conformation, and
the evaluator classifies the fitness of the proposed structure based on a
given ruleset. Normally, this is the energy of the system according to a
given forcefield. The proposition is then accepted or rejected based on a
criterium (normally the Metropolis criterium). In this criteria, structures
of lower energy are automatically accepted, where higher energy structures
may or not be accepted depending on the energy difference between the previous
solution and the energy of the system (higher energies give higher probabilities
of worst solutions being accepted), thus allowing the simulation to overcome
small bumps in the free-energy landscape. The system is allowed to roam the local
conformational neighbourhood for a finite number of steps, visiting favorable
structures more frequently.

# Convergence criteria
1. Number of `n_steps` performed;
---------------------------------------------------------------------------- =#

@doc raw"""
    DriverConfig(sampler::Abstract.Sampler, evaluator::Abstract.Evaluator, [, temperature::Union{Function, Float64} = 0.0, n_steps::Int64 = 0, callbacks::Vector{Abstract.CallbackObject} = []])

Define the runtime parameters for the Monte Carlo simulation.
No `sampler` movement is performed by default, since n_steps = 0.

# Arguments
- `sampler::Abstract.Sampler`: Responsible for generating a new structure to be evaluated. See [`Sampler`](@ref Mutators).
- `evaluator::Abstract.Evaluator`: Responsible for evaluating the system energy. See [`Evaluator`](@ref Forcefield).
- `temperature::Union{Function, Float64}`: (Optional) Temperature of the system, determines acceptance in the Metropolis algorithm.
This parameter is always a function who should follow the following signature:
```
function (n::Int64)::Float64
```
When a Float64 is passed, by default, ProtoSyn will construct an unnamed function who always returns that Float64 number every step (Default: 0.0)
- `n_steps::Int64`: (Optional) Total amount of steps to be performed (Default: 0).
- `callbacks::Vector{Abstract.CallbackObject}`: (Optional) Vector of [`CallbackObject`](@ref Common)s.

# Examples
```julia-repl
julia> Drivers.MonteCarlo.DriverConfig(sampler = my_sampler!, evaluator = my_evaluator!, temperature = 10.0, n_steps = 1000)
MonteCarlo.Driver(sampler=my_sampler!, evaluator=my_evaluator!, temperature=10.0, n_steps=1000, callbacks=[])

julia> Drivers.MonteCarlo.DriverConfig(my_sampler!, my_evaluator!)
MonteCarlo.Driver(sampler=my_sampler!, evaluator=my_evaluator!, temperature=0.0, n_steps=0, callbacks=[])
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
    callbacks::Vector{T}                     # Default: empty

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


@doc raw"""
    DriverState(;step::Int64 = 0, ac_count::Int64 = -1, temperature::Float64 = -1.0, completed::Bool = false)

Define the runtime parameters given by the Monte Carlo simulation.

# Arguments
- `step::Int64`: Current step of the simulation.
- `ac_count::Int64`: Current acceptance count of the simulation. Gives acceptance ration when divided by `step`.
- `temperature::Float64`: Current temperature of the system.
- `completed::Bool`: Current completness status of the system.

# Examples
```julia-repl
julia> Drivers.MonteCarlo.DriverState()
ProtoSyn.Drivers.MonteCarlo.DriverState
   step = 0
   ac_count = -1
   temperature = -1.0
   completed = false
```
See also: [`run!`](@ref)
"""
Base.@kwdef mutable struct DriverState <: Abstract.DriverState
    
    # Parameter:         # Default:
    step::Int64          = 0
    ac_count::Int64      = -1
    temperature::Float64 = -1.0
    completed::Bool      = false
end


# ----------------------------------------------------------------------------------------------------------
#                                                   RUN

@doc raw"""
    run!(state::Common.State, driver_config::DriverConfig)

Run the main body of the driver. Create a new conformation based on `driver_config.sampler`, evaluate the new conformation energy using `driver_config.evaluator`,
accepting it or not depending on the `driver_state.temperature` in a Metropolis algorithm (updated based on `driver_config.anneal_fcn`). This Monte Carlo process 
is repeated for `driver_config.n_steps`, saving the accepted structures to `state` and calling all the `driver_config.callbacks`. 

# Arguments
- `state::Common.State`: Current state of the system to be modified.
- `driver_config::DriverConfig`: Defines the parameters for the MonteCarlo simulation. See [`DriverConfig`](@ref).

# Examples
```julia-repl
julia> Drivers.MonteCarlo.run!(state, mc_driver)
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
    
    # Instantiate a new DriverState object.
    # By default, no optimization step has yet been taken
    # apart from calculating the energy and forces for the input state
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
        end # end if
        
        # update driver state and call calback functions (if any)
        Common.@cbcall driver_config.callbacks state driver_state
        if driver_config.sampler.tune! != nothing
            driver_config.sampler.tune!(driver_config.sampler.mutators, driver_state)
        end # end if

    end # end while
    #endregion

    # update driver state and return
    driver_state.completed = true

    return  driver_state
    
end #end function

end #end module