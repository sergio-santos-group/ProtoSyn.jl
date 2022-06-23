using ProtoSyn.Calculators: EnergyFunction
using ProtoSyn.Mutators: AbstractMutator
using Printf

"""
    MonteCarloState{T <: AbstractFloat}(step::Int = 0, converged::Bool = false, completed::Bool = false, stalled::Bool = false, acceptance_count = 0, temperature::T = T(0.0))
    MonteCarloState{T <: AbstractFloat}()

A [`MonteCarloState`](@ref) instance is a `DriverState` that keeps track of a
[`MonteCarlo`](@ref) simulation status. Besides the default `:step`,
`:converged`, `:completed` and `:stalled` fields, a [`MonteCarloState`](@ref)
instance adds the following tracking fields specific for [`MonteCarlo`](@ref)
simulations:

* `acceptance_count::Int` - The number of accepted [`Pose`](@ref) changes by the Metropolis criterion;
* `temperature::T` - The current temperature of the simulation.

# See also
[`MonteCarlo`](@ref)

# Examples
```jldoctest
julia> ProtoSyn.Drivers.MonteCarloState{Float64}()
ProtoSyn.Drivers.MonteCarloState{Float64}(0, false, false, false, 0, 0.0)
```
"""
Base.@kwdef mutable struct MonteCarloState{T <: AbstractFloat} <: DriverState
    step::Int             = 0
    converged::Bool       = false
    completed::Bool       = false
    stalled::Bool         = false
    acceptance_count::Int = 0
    temperature::T        = T(0.0)
end


"""
    MonteCarlo(eval!::Union{Function, EnergyFunction}, sample!::Union{Function, AbstractMutator, Driver}, callback::Opt{Callback}, max_steps::Int, temperature::Function)

A [`MonteCarlo`](@ref) `Driver` instance. As such, this object is callable as a
_functor_ with the following signature:

```julia
(driver::MonteCarlo)(pose::Pose)
```

This simulation performs a Markov Chain Monte Carlo algorithm (MCMC), more
especifically, the Metropolis-Hastings method. Starting from an initial
[`Pose`](@ref) `pose` conformation, new changes are imposed by calling the
given `sample!` method (this is an `AbstractMutator`, `Driver` or a custom
sampling function). Note that the [`Pose`](@ref) is synched (using the
[`sync!`](@ref) method) after each `sample!` call. The two [`Pose`](@ref)
instances are evaluated by calling the `eval!` method (an
[`EnergyFunction`](@ref) or custom evaluator function). Following the
[Metropolis Criterion](https://en.wikipedia.org/wiki/Metropolis%E2%80%93Hastings_algorithm),
the new [`Pose`](@ref) is accepted if the evaluated energy is lower than the
previous [`Pose`](@ref), or, otherwise, with a given probability based on the
`temperature` of the simulation:

\$\\;\\;\\;\\;\\;\\; rand()<e^{-\\frac{E_{n+1}-E_{n}}{T}} \\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\\;\$

Where \$E_{n+1}\$ is the energy of new [`Pose`](@ref), \$E_{n}\$ is the current
[`Pose`](@ref) energy and \$T\$ is current `temperature` of the simulation. As
such, higher `temperature` values correlate with a higher chance to accept
higher energy [`Pose`](@ref) instances (in order to overcome energy barriers).
If rejected, the previous [`Pose`](@ref) is recovered for a new attempt. This
simulation loop (sample! -> eval! -> metropolis criterion) continues for a total
of `max_steps`. Each step, a call to an optional [`Callback`](@ref) `callback`
is performed. A companion [`MonteCarloState`](@ref) `DriverState` instance is
also updated each step and provided to the [`Callback`](@ref) `callback`.

# Fields
* `eval!::Union{Function, EnergyFunction}` - The evaluator [`EnergyFunction`](@ref) or custom function, receives a [`Pose`](@ref) `pose` as the single argument;
* `sample!::Union{Function, AbstractMutator, Driver}` - The sampler method, receives a [`Pose`](@ref) `pose` as the single input argument;
* `callback::Opt{Callback}` - An optional [`Callback`](@ref) instance, receives two input arguments: the current [`Pose`](@ref) `pose` and the current `DriverState` `driver_state`;
* `max_steps::Int` - The total number of simulation steps to be performed;
* `temperature::Function` - The temperature function thermostat (see [Available thermostats](@ref)), receives the current step (as an `Int`) as the single input argument.

# See also
[`MonteCarloState`](@ref)

# Examples
```jldoctest
julia> monte_carlo = ProtoSyn.Drivers.MonteCarlo(
                  energy_function,
                  dihedral_mutator,
                  ProtoSyn.Common.default_energy_step_callback(1),
                  10,
                  ProtoSyn.Drivers.get_linear_quench(1.0, 10))
⚒  Monte Carlo Driver:
 ├──  ●  Evaluator:
 |    └── 🗲  Energy Function (4 components):
 |        +----------------------------------------------------------------------+
 |        | Index | Component name                                | Weight (α)   |
 |        +----------------------------------------------------------------------+
 |        | 1     | TorchANI_ML_Model                             |      1.000   |
 |        | 2     | Caterpillar_Solvation                         |      0.010   |
 |        | 3     | Bond_Distance_Restraint                       |      1.000   |
 |        | 4     | Cα-Cα_Clash_Restraint                         |    100.000   |
 |        +----------------------------------------------------------------------+
 |   
 ├──  ● Sampler:
 |    └── ⚯  Dihedral Mutator:
 |        +----------------------------------------------------------------------+
 |        | Index | Field                       | Value                          |
 |        +----------------------------------------------------------------------+
 |        | 1     | angle_sampler               | Function randn                 |
 |        | 2     | p_mut                       | 0.0100                         |
 |        | 3     | step_size                   | 0.5000                         |
 |        +----------------------------------------------------------------------+
 |         ● Selection: Set
 |         └── FieldSelection › Atom.name = r"C|N"
 |   
 ├──  ● Callback:
 |    └── ✉  Callback:
 |        +----------------------------------------------------------------------+
 |        | Index | Field                     | Value                            |
 |        +----------------------------------------------------------------------+
 |        | 1     | Event                     | energy_step                      |
 |        | 2     | Frequency                 | 1                                |
 |        +----------------------------------------------------------------------+
 |   
 └──  ● Settings:
       Max steps: 10
     Temperature: linear_quench
```
"""
mutable struct MonteCarlo <: Driver
    eval!::Union{Function, EnergyFunction}
    sample!::Union{Function, AbstractMutator, Driver}
    callback::Opt{Callback}
    max_steps::Int
    temperature::Function
end


function (driver::MonteCarlo)(pose::Pose)
    
    T = eltype(pose.state)
    driver_state = MonteCarloState{T}()
    driver_state.temperature = driver.temperature(0)

    previous_energy = driver.eval!(pose)
    previous_state  = copy(pose)
    driver.callback !== nothing && driver.callback(pose, driver_state)
    
    while driver_state.step < driver.max_steps
        
        driver.sample!(pose)
        sync!(pose)
        energy = driver.eval!(pose)
        
        n = rand()
        driver_state.temperature = driver.temperature(driver_state.step)
        m = exp((-(energy - previous_energy)) / driver_state.temperature)
        if (energy < previous_energy) || (n < m)
            # println("Accepted")
            previous_energy = energy
            ProtoSyn.recoverfrom!(previous_state, pose) # If copy, the chain is broken
            driver_state.acceptance_count += 1
        else
            # println("Rejected")
            ProtoSyn.recoverfrom!(pose, previous_state) # If copy, the chain is broken
        end

        driver_state.step += 1
        driver.callback !== nothing && driver.callback(pose, driver_state)
    end

    driver_state.completed = true
    return pose
end

# * Show -----------------------------------------------------------------------

function Base.show(io::IO, mcd::MonteCarlo, level_code::Opt{LevelCode} = nothing)
    init_level_code = level_code === nothing ? LevelCode() : level_code
    init_lead       = ProtoSyn.get_lead(level_code)
    init_inner_lead = ProtoSyn.get_inner_lead(level_code)

    println(io, init_lead*"⚒  Monte Carlo Driver:")

    level_code = vcat(init_level_code, 3)
    lead       = ProtoSyn.get_lead(level_code)
    inner_lead = ProtoSyn.get_inner_lead(level_code)

    if typeof(mcd.eval!) === ProtoSyn.Calculators.EnergyFunction
        println(io, lead*" ●  Evaluator:")
        Base.show(io, mcd.eval!, vcat(level_code, 4))
    else
        println(io, lead*" ●  Evaluator: $(mcd.eval!)")
    end

    println(io, init_inner_lead*level_code.code_table[1])

    if isa(mcd.sample!, AbstractMutator)
        println(io, lead*" ● Sampler:")
        Base.show(io, mcd.sample!, vcat(level_code, 4))
    elseif isa(mcd.sample!, Driver)
        println(io, lead*" ● Sampler:")
        Base.show(io, mcd.sample!, vcat(level_code, 4))
    else
        println(io, lead*" ●  Sampler: $(mcd.sample!)")
    end

    println(io, init_inner_lead*level_code.code_table[1])

    if mcd.callback !== nothing
        println(io, lead*" ● Callback:")
        Base.show(io, mcd.callback, vcat(level_code, 4))
    else
        println(io, lead*" ○  Callback: Not set")
    end

    level_code = vcat(init_level_code, 4)
    lead       = ProtoSyn.get_lead(level_code)
    inner_lead = ProtoSyn.get_inner_lead(level_code)

    println(io, init_inner_lead*level_code.code_table[1])

    println(io, lead*" ● Settings:")
    println(io, inner_lead*"  Max steps: $(mcd.max_steps)")
    println(io, inner_lead*"Temperature: $(mcd.temperature)")
end