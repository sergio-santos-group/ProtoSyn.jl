using ProtoSyn.Calculators: EnergyFunction
using ProtoSyn.Mutators: AbstractMutator

"""
    ILSState{T <: AbstractFloat}(step::Int = 0, converged::Bool = false, completed::Bool = false, stalled::Bool = false, acceptance_count = 0, temperature::T = T(0.0))
    ILSState{T <: AbstractFloat}()

A [`ILSState`](@ref) instance is a `DriverState` that keeps track of a
[`ILS`](@ref) simulation status. Besides the default `:step`,
`:converged`, `:completed` and `:stalled` fields, a [`ILSState`](@ref)
instance adds the following tracking fields specific for [`ILS`](@ref)
simulations:

* `acceptance_count::Int` - The number of accepted [`Pose`](@ref) changes by the Metropolis criterion;
* `temperature::T` - The current temperature of the simulation.

# See also
[`ILS`](@ref)

# Examples
```jldoctest
julia> ProtoSyn.Drivers.ILSState{Float64}()
ProtoSyn.Drivers.ILSState{Float64}(0, false, false, false, 0, 0.0)
```
"""
Base.@kwdef mutable struct ILSState{T <: AbstractFloat} <: DriverState
    step::Int        = 0
    converged::Bool  = false
    completed::Bool  = false
    stalled::Bool    = false
    acceptance_count = 0
    temperature::T   = T(0.0)
end


"""
    ILS(eval!::Union{Function, EnergyFunction}, jump!::Union{Function, AbstractMutator, Driver}, inner_driver!::Driver, callback::Opt{Callback}, max_steps::Int, temperature::Function)

An [`ILS`](@ref) `Driver` instance. As such, this object is callable as a
_functor_ with the following signature:

```julia
(driver::ILS)(pose::Pose)
```

This simulation performs an
[Iterated Local Search](https://en.wikipedia.org/wiki/Iterated_local_search)
algorithm. In this type of simulation, an `inner_driver!` takes the system from
a state to another, possibly optimizing the structure into a local minimum. As
an example, this can be a [`MonteCarlo`](@ref) or a [`SteepestDescent`](@ref)
`Driver`. The resulting structure is evaluated (according to an `eval!`
[`EnergyFunction`](@ref) instance or custom function) and accepted or rejected
according to the Metropolis Criterion, based on the current `temperature` of the
simulation (see [`MonteCarlo`](@ref) for more information regarding the
Metropolis Criterion). In either case, the accepted structure or the recovered
previous state structure is subjected to a `jump!`, that is, a large
conformational change, taking the system to a new conformational space. This
simulation loop runs for `max_steps` iteration. Each step, a call to an optional
[`Callback`](@ref) `callback` is performed. A companion [`ILSState`](@ref)
`DriverState` instance is also updated each step and provided to the
[`Callback`](@ref) `callback`.

# Fields
* `eval!::Union{Function, EnergyFunction}` - The evaluator [`EnergyFunction`](@ref) or custom function, receives a [`Pose`](@ref) `pose` as the single argument;
* `jump!::Union{Function, AbstractMutator, Driver}` - The jump method, receives a [`Pose`](@ref) `pose` as the single input argument, should introduce a relative high change in the conformation;
* `inner_driver!::Driver` - The inner loop `Driver`;
* `callback::Opt{Callback}` - An optional [`Callback`](@ref) instance for the outer loop, receives two input arguments: the current [`Pose`](@ref) `pose` and the current `DriverState` `driver_state`;
* `max_steps::Int` - The total number of simulation steps to be performed in the outer loop;
* `temperature::Function` - The temperature function thermostat for the outer loop (see [Available thermostats](@ref)), receives the current step (as an `Int`) as the single input argument.

# See also
[`ILSState`](@ref)

# Examples
```jldoctest
julia> ils = ProtoSyn.Drivers.ILS(energy_function, dihedral_mutator, monte_carlo, cb, 1000, ProtoSyn.Drivers.get_linear_quench(1.0, 1000))
⚒  ILS Driver:
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
 ├──  ● Jump:
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
 ├──  ● Inner Driver:
 |    └── ⚒  Monte Carlo Driver:
 |         ├──  ●  Evaluator:
 |         |    └── 🗲  Energy Function (4 components):
 |         |        +----------------------------------------------------------------------+
 |         |        | Index | Component name                                | Weight (α)   |
 |         |        +----------------------------------------------------------------------+
 |         |        | 1     | TorchANI_ML_Model                             |      1.000   |
 |         |        | 2     | Caterpillar_Solvation                         |      0.010   |
 |         |        | 3     | Bond_Distance_Restraint                       |      1.000   |
 |         |        | 4     | Cα-Cα_Clash_Restraint                         |    100.000   |
 |         |        +----------------------------------------------------------------------+
 |         |   
 |         ├──  ● Sampler:
 |         |    └── ⚯  Dihedral Mutator:
 |         |        +----------------------------------------------------------------------+
 |         |        | Index | Field                       | Value                          |
 |         |        +----------------------------------------------------------------------+
 |         |        | 1     | angle_sampler               | Function randn                 |
 |         |        | 2     | p_mut                       | 0.0100                         |
 |         |        | 3     | step_size                   | 0.5000                         |
 |         |        +----------------------------------------------------------------------+
 |         |         ● Selection: Set
 |         |         └── FieldSelection › Atom.name = r"C|N"
 |         |   
 |         ├──  ○  Callback: Not set
 |         |   
 |         └──  ● Settings:
 |               Max steps: 10
 |             Temperature: linear_quench
 |   
 ├──  ● Callback:
 |    └── ✉  Callback:
 |        +----------------------------------------------------------------------+
 |        | Index | Field                     | Value                            |
 |        +----------------------------------------------------------------------+
 |        | 1     | Event                     | energy_step_frame                |
 |        | 2     | Frequency                 | 10                               |
 |        +----------------------------------------------------------------------+
 |   
 └──  ● Settings:
       Max steps: 1000
     Temperature: linear_quench
```
"""
mutable struct ILS <: Driver
    eval!::Union{Function, EnergyFunction}
    jump!::Union{Function, AbstractMutator, Driver}
    inner_driver!::Driver
    callback::Opt{Callback}
    max_steps::Int
    temperature::Function # Takes step, returns temperature at step
end


function (driver::ILS)(pose::Pose)
    
    T = eltype(pose.state)
    driver_state = ILSState{T}()
    driver_state.temperature = driver.temperature(0)
    
    previous_state  = copy(pose)
    previous_energy = driver.eval!(pose)
    driver.callback !== nothing && driver.callback(pose, driver_state)
    
    while driver_state.step < driver.max_steps
            
        driver.inner_driver!(pose)
        sync!(pose)
        energy = driver.eval!(pose)
        
        n = rand()
        driver_state.temperature = driver.temperature(driver_state.step)
        m = exp((-(energy - previous_energy)) / driver_state.temperature)
        if (energy < previous_energy) || (n < m)
            previous_energy = energy
            previous_state  = copy(pose)
            driver_state.acceptance_count += 1
        else
            ProtoSyn.recoverfrom!(pose, previous_state)
        end

        driver_state.step += 1
        driver.callback !== nothing && driver.callback(pose, driver_state)

        # Jump
        if driver_state.step < driver.max_steps
            driver.jump!(pose)
            sync!(pose)
        end
    end

    driver_state.completed = true
    return pose
end

# * Show -----------------------------------------------------------------------

function Base.show(io::IO, ils::ILS, level_code::Opt{LevelCode} = nothing)
    init_level_code = level_code === nothing ? LevelCode() : level_code
    init_lead       = ProtoSyn.get_lead(level_code)
    init_inner_lead = ProtoSyn.get_inner_lead(level_code)

    println(io, init_lead*"⚒  ILS Driver:")

    level_code = vcat(init_level_code, 3)
    lead       = ProtoSyn.get_lead(level_code)
    inner_lead = ProtoSyn.get_inner_lead(level_code)

    if typeof(ils.eval!) === ProtoSyn.Calculators.EnergyFunction
        println(io, lead*" ●  Evaluator:")
        Base.show(io, ils.eval!, vcat(level_code, 4))
    else
        println(io, lead*" ●  Evaluator: $(ils.eval!)")
    end

    println(io, init_inner_lead*level_code.code_table[1])

    if isa(ils.jump!, AbstractMutator)
        println(io, lead*" ● Jump:")
        Base.show(io, ils.jump!, vcat(level_code, 4))
    elseif isa(ils.jump!, Driver)
        println(io, lead*" ● Jump:")
        Base.show(io, ils.jump!, vcat(level_code, 4))
    else
        println(io, lead*" ●  Jump: $(ils.jump!)")
    end

    println(io, init_inner_lead*level_code.code_table[1])

    println(io, lead*" ● Inner Driver:")
    Base.show(io, ils.inner_driver!, vcat(level_code, 4))

    println(io, init_inner_lead*level_code.code_table[1])

    if ils.callback !== nothing
        println(io, lead*" ● Callback:")
        Base.show(io, ils.callback, vcat(level_code, 4))
    else
        println(io, lead*" ○  Callback: Not set")
    end

    level_code = vcat(init_level_code, 4)
    lead       = ProtoSyn.get_lead(level_code)
    inner_lead = ProtoSyn.get_inner_lead(level_code)

    println(io, init_inner_lead*level_code.code_table[1])

    println(io, lead*" ● Settings:")
    println(io, inner_lead*"  Max steps: $(ils.max_steps)")
    println(io, inner_lead*"Temperature: $(ils.temperature)")
end