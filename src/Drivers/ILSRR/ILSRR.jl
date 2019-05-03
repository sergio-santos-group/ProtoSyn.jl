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
Base.@kwdef mutable struct DriverConfig{F <: Function, G <: Function, H <: Function} <: Drivers.AbstractDriverConfig

    inner_driver::Drivers.AbstractDriverConfig
    evaluator!::F
    perturbator!::G
    anneal_fcn::H
    n_steps::Int64 = 0
    continue_after_n_attemps::Int64 = 0

end
function DriverConfig(inner_driver::Drivers.AbstractDriverConfig, evaluator!::F, perturbator!::G, temperature::Float64 = 0.0) where {F <: Function, G <: Function}
    return DriverConfig(inner_driver = inner_driver, evaluator! = evaluator!, perturbator! = perturbator!, anneal_fcn = (n::Int64)->temperature)
end

# TODO: Documentation
Base.@kwdef mutable struct DriverState <: Drivers.AbstractDriverState
    step::Int            = 0
    fj_count::Int        = 0 # Failed jumps count
    temperature::Float64 = -1.0
    completed::Bool      = false
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

    save_inner_best = @Common.callback 1 function _save_inner_best(st::Common.State, dr_state::Drivers.AbstractDriverState, dr_config::Drivers.AbstractDriverConfig)
        if st.energy.eTotal < inner_best.energy.eTotal
            copy!(inner_best.xyz, st.xyz)
            copy!(inner_best.energy, st.energy)
        end
    end

    function set_homebase()
        # printstyled(@sprintf("(ILSRR) New homebase defined: ⚡E: %10.3e (old) ▶️ %10.3e (new)\n", homebase.energy.eTotal, state.energy.eTotal), color = :green)
        copy!(homebase.xyz, state.xyz)
        copy!(homebase.energy, state.energy)
        driver_state.fj_count = 0
    end

    driver_state = DriverState()
    driver_config.evaluator!(state, false)

    inner_best   = Common.State(state.size)
    copy!(inner_best.xyz, state.xyz)
    copy!(inner_best.energy, state.energy)

    homebase     = Common.State(state.size)
    copy!(homebase.xyz, state.xyz)
    copy!(homebase.energy, state.energy)

    initial_state = Common.State(state.size)
    copy!(initial_state.xyz, state.xyz)
    copy!(initial_state.energy, state.energy)
    
    best_energy   = Inf 

    while driver_state.step < driver_config.n_steps
        # println(@sprintf("\n(%5s) %12s \n%s", "ILSRR", @sprintf("Step: %4d", step), "-"^150))

        copy!(inner_best.xyz, st.xyz)
        copy!(inner_best.energy, st.energy)
        driver_config.inner_cycle_driver.run!(state, driver_config.inner_driver, save_inner_best) #!!!!!!
        state = deepcopy(inner_best)
        
        @Common.cbcall driver.callbacks..., callbacks... step state driver
        
        if state.energy.eTotal < best_energy
            best_energy = state.energy.eTotal
            set_homebase()
        else
            if state.energy.eTotal < homebase.energy.eTotal || (rand() < exp(-(state.energy.eTotal - homebase.energy.eTotal) / driver.temperature)) # Metropolis Criteria
                set_homebase()
            else
                failed_jumps_count += 1
                if driver.continue_after_n_attemps > 0 && failed_jumps_count >= driver.continue_after_n_attemps
                    printstyled(@sprintf("(ILSRR) Reseting because %2d consecutive jumps failed to produce a new homebase\n", failed_jumps_count), color = :red)
                    state    = deepcopy(initial_state)
                    homebase = deepcopy(initial_state)
                    failed_jumps_count = 0
                    continue
                end
                printstyled(@sprintf("(ILSRR) Recovering to previous homebase (x%2d): ⚡E: %10.3e (actual) ▶️ %10.3e (new)\n", failed_jumps_count, state.energy.eTotal, homebase.energy.eTotal), color = 9)
                state = deepcopy(homebase)
            end
            if step != driver.n_steps
                driver.perturbator!(state)
                driver.evaluator!(state, false)
            end
        end
    end
end

end