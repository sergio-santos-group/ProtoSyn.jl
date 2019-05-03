module ILSRR

using ..Aux
using ..Common
using ..Drivers
using Printf


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
mutable struct DriverConfig{F <: Function, G <: Function}
    #TO DO: Documentation
    n_steps::Int
    temperature::Float64
    perturbator!::G
    inner_driver!::F
    inner_driver_config::Union{Drivers.AbstractDriverConfig,Nothing}
end

# Base.@kwdef mutable struct DriverConfig{F <: Function, G <: Function, H <: Function} <: Drivers.AbstractDriverConfig

#     inner_driver::Drivers.AbstractDriverConfig
#     evaluator!::F
#     perturbator!::G
#     anneal_fcn::H
#     n_steps::Int64 = 0
#     continue_after_n_attemps::Int64 = 0

# end
# function DriverConfig(inner_driver::Drivers.AbstractDriverConfig, evaluator!::F, perturbator!::G, temperature::Float64 = 0.0) where {F <: Function, G <: Function}
#     return DriverConfig(inner_driver = inner_driver, evaluator! = evaluator!, perturbator! = perturbator!, anneal_fcn = (n::Int64)->temperature)
# end
# Base.show(io::IO, b::DriverConfig) = print(io, "ILSRR.DriverConfig(inner_cycle_driver=$(b.inner_cycle_driver), evaluator=$(string(b.evaluator!)), perturbator=$(string(b.perturbator!)), anneal_fcn=$(string(b.anneal_fcn)), n_steps=$(b.n_steps), continue_after_n_attemps=$(b.continue_after_n_attemps))")


#TO DO: Documentation
Base.@kwdef mutable struct DriverState <: Drivers.AbstractDriverState
    step::Int64        = 0
    n_stalls::Int      = 0
    completed::Bool    = false
    stalled::Bool      = false
    best_state::Common.State
    home_state::Common.State
end


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
function run!(state::Common.State, driver_config::DriverConfig, callbacks::Common.CallbackObject...)

    # save_inner_best = @Common.callback 1 function cb_save(step::Int64, st::Common.State, dr::Drivers.AbstractDriverConfig, args...)
    #     if st.energy.eTotal < inner_best.energy.eTotal
    #         inner_best = deepcopy(st)
    #     end
    # end

    # save_inner_best = Drivers.@callback 1 function(st)
    #     if st.energy.total < inner_best.energy.total
    #         # just copy coordinates and energy:
    #         #   inner_best.energy = st.energy
    #         #   inner_best.xyz = st.xyz
    #         Common.@copy inner_best st energy xyz
    #     end
    # end

    # function set_homebase()
    #     printstyled(@sprintf("(ILSRR) New homebase defined: ⚡E: %10.3e (old) ▶️ %10.3e (new)\n", homebase.energy.eTotal, state.energy.eTotal), color = :green)
    #     Common.@copy homebase state energy xyz
    #     # copy!(homebase, state)
    #     failed_jumps_count = 0
    # end

    # step::Int64 = 1
    # failed_jumps_count::Int64 = 0
    # driver.evaluator!(state, false)
    # inner_best    = deepcopy(state)
    # homebase      = deepcopy(state)
    # initial_state = deepcopy(state)
    # best_energy   = Inf 

    # for step in 1:driver.n_steps
    #     println(@sprintf("\n(%5s) %12s \n%s", "ILSRR", @sprintf("Step: %4d", step), "-"^150))

    #     inner_best = deepcopy(initial_state)
    #     driver.inner_cycle_driver.run!(state, driver.inner_cycle_driver, save_inner_best)
    #     state = deepcopy(inner_best)
        
    #     @Common.cbcall driver.callbacks..., callbacks... step state driver
        
    #     if state.energy.total < best_energy
    #         best_energy = state.energy.total
    #         set_homebase()
    #     else
    #         ΔE = state.energy.total - homebase.energy.total
    #         if (ΔE < 0.0) || (rand() < exp(-ΔE/driver.temperature)) # Metropolis Criteria
    #             set_homebase()
    #         else
    #             failed_jumps_count += 1
    #             if driver.continue_after_n_attemps > 0 && failed_jumps_count >= driver.continue_after_n_attemps
    #                 printstyled(@sprintf("(ILSRR) Reseting because %2d consecutive jumps failed to produce a new homebase\n", failed_jumps_count), color = :red)
    #                 state    = deepcopy(initial_state)
    #                 homebase = deepcopy(initial_state)
    #                 failed_jumps_count = 0
    #                 continue
    #             end
    #             printstyled(@sprintf("(ILSRR) Recovering to previous homebase (x%2d): ⚡E: %10.3e (actual) ▶️ %10.3e (new)\n", failed_jumps_count, state.energy.eTotal, homebase.energy.eTotal), color = 9)
    #             state = deepcopy(homebase)
    #         end
    #         if step != driver.n_steps
    #             driver.perturbator!(state)
    #             driver.evaluator!(state, false)
    #         end
    #     end
    # end

    #--------------------------------------------------------------------------
    n_stalls = 0
    driver_state = DriverState()
    

    inner_driver! = driver_config.inner_driver
    inner_driver_config = driver_config.inner_driver_config
    
    let n_steps=inner_driver_config.n_steps
        inner_driver_config.n_steps = 0
        inner_driver!(state, inner_driver_config)
        inner_driver_config.n_steps = n_steps
    end

    driver_state.best_state = Common.State(state)
    driver_state.home_state = Common.State(state)
    driver_state.completed = driver_state.step == driver_config.n_steps

    @Common.cbcall callbacks state driver_state driver_config
    

    #region MAINLOOP
    while !(driver_state.completed || driver_state.stalled)
        
        # this driver should make multiple small tweaks
        # to the state
        inner_driver!(state, inner_driver_config, callbacks)
        
        driver_state.step += 1
        @Common.cbcall callbacks state driver_state driver_config

        if state.energy.total < best.energy.total
            # save this state as the best and make
            # it the new homebase
            Common.@copy best state energy xyz
            Common.@copy homebase state energy xyz
            n_stalls = 0
        else
            # otherwise, a new homebase may be created according
            # to the Metropolis criterium
            ΔE = state.energy.total - homebase.energy.total
            if (ΔE <= 0.0) || (rand() < exp(-ΔE/driver.temperature))
                Common.@copy homebase state energy xyz
                n_stalls = 0
            else
                # if the criterium was not accepted, revert
                # to the homebase
                Common.@copy state homebase energy xyz
                n_stalls += 1
            end
        end

        driver_state.stalled = n_stalls == driver_config.stall_limit
        driver_state.completed = driver_state.step == driver_config.n_steps
        
        # make a large perturbation to the state
        if !(driver_state.completed || driver_state.stalled)
            driver_config.perturbator!(state)
        end

    end
    #endregion

    # before returning, save the best state
    if best.energy.total < state.energy.total
        Common.@copy state best energy xyz
    end

    return driver_state

end # end function


end # end module