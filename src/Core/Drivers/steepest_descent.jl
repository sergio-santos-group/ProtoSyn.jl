using ProtoSyn
using ProtoSyn.Calculators: EnergyFunction

"""
    SteepestDescentState{T <: AbstractFloat}(step::Int = 0, converged::Bool = false, completed::Bool = false, stalled::Bool = false, stepsize::T = T(1), max_force::Tuple{T,Int} = (T(-1),0))
    SteepestDescentState{T <: AbstractFloat}()

A [`SteepestDescentState`](@ref) instance is a `DriverState` that keeps track of
a [`SteepestDescent`](@ref) simulation status. Besides the default `:step`,
`:converged`, `:completed` and `:stalled` fields, a
[`SteepestDescentState`](@ref) instance adds the following tracking fields
specific for [`SteepestDescent`](@ref) simulations:

* `stepsize::T` - The current learning rate of the simulation;
* `max_force::Tuple{T,Int}` - The maximum force felt by the system and the atom index where the force is felt.

# See also
[`SteepestDescent`](@ref)

# Examples
```jldoctest
julia> ProtoSyn.Drivers.SteepestDescentState{Float64}()
ProtoSyn.Drivers.SteepestDescentState{Float64}(0, false, false, false, 1.0, (-1.0, 0))
```
"""
Base.@kwdef mutable struct SteepestDescentState{T <: AbstractFloat} <: DriverState
    step::Int          = 0
    converged::Bool    = false
    completed::Bool    = false
    stalled::Bool      = false

    stepsize::T        =  T(1)
    max_force::Tuple{T,Int} = (T(-1),0)
end


"""
    SteepestDescent(eval!::Union{Function, EnergyFunction}, callback::Opt{Callback}, max_steps::Int, force_tolerance::Float64, max_displacement::Float64)

A [`SteepestDescent`](@ref) `Driver` instance. As such, this object is callable
as a _functor_ with the following signature:

```julia
(driver::SteepestDescent)(pose::Pose)
```

This simulation employs the calculation of forces by an evaluating function or
[`EnergyFunction`](@ref) instance (`eval!`) to relax and optimize a given
[`Pose`](@ref) `pose`. This process will run for a maximum of `max_steps`
iterations, or until the simulation has stalled or converged. Convergence is
achieved when the maximum force felt on the system is bellow a given threshold
`force_tolerance`, while a simulation is stalled when the current scaling factor
is bellow the machine precision. In each step, the displacement is calculated as
a factor of the force applied to an atom multiplied by a step size (or learning
rate), up to a maximum `max_displacement`. As an heuristic for the definition of
this step size, this value increases while when the system's energy is lowering,
and rapidly decreases if an higher energy [`State`](@ref) is generated. This
allows the simulation to traverse for some degree of energy barriers. Each step,
a call to an optional [`Callback`](@ref) `callback` is performed. A companion
[`SteepestDescentState`](@ref) `DriverState` instance is also updated each step
and provided to the [`Callback`](@ref) `callback`.

# Fields
* `eval!::Union{Function, EnergyFunction}` - The evaluator [`EnergyFunction`](@ref) or custom function, receives a [`Pose`](@ref) `pose` as the single argument;
* `callback::Opt{Callback}` - An optional [`Callback`](@ref) instance, receives two input arguments: the current [`Pose`](@ref) `pose` and the current `DriverState` `driver_state`;
* `max_steps::Int` - The total number of simulation steps to be performed;
* `force_tolerance::Float64` - The minimum force tolerated by the simulation, any [`State`](@ref) with a max force below this threshold will converge the simulation;
* `max_displacement::Float64` - The maximum displacement each step (a.u.).

# See also
[`SteepestDescentState`](@ref)

# Examples
```jldoctest
julia> sd = ProtoSyn.Drivers.SteepestDescent(energy_function, cb, 1000, 0.001, 0.1)
⚒  Steepest Descent Driver:
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
      Force tolerance: 0.001
     Max displacement: 0.1
```
"""
mutable struct SteepestDescent <: Driver
    eval!::Union{Function, EnergyFunction}
    callback::Opt{Callback}
    max_steps::Int
    force_tolerance::Float64
    max_displacement::Float64
end

SteepestDescent(::Type{T}) where {T <: AbstractFloat} = begin
    zero_eval = (pose::Pose, update_forces::Bool = False) -> return 0.0
    return SteepestDescent(zero_eval, nothing, 0, T(1), T(0.1))
end

SteepestDescent() = begin
    return SteepestDescent(ProtoSyn.Units.defaultFloat)
end


function (driver::SteepestDescent)(pose::Pose)
  
    # SteepestDescent requires a synched pose
    sync!(pose)

    T = eltype(pose.state)
    energy = driver.eval!(pose, update_forces_overwrite = true)

    # Instantiate a new DriverState object. By default, no optimization step has
    # yet been taken apart from calculating the energy and forces for the input
    # state
    driver_state = SteepestDescentState{T}()
    driver_state.max_force = ProtoSyn.atmax(pose.state, :f)
    
    # If max force is already below the requested force tolerance, simply return
    # the current driver state
    driver_state.converged = driver_state.max_force[1] < driver.force_tolerance

    if driver_state.converged
        println("Steepest descent converged to minimum force tolerance ($(driver_state.max_force[1]) < $(driver.force_tolerance)) at step $(driver_state.step).")
        return driver_state
    end
    
    # Initial callbacks
    driver.callback !== nothing && driver.callback(pose, driver_state)

    # Initialize γ and backup copy
    γ    = T(1.0)
    γinc = T(1.05)
    γdec = T(0.50) 
    γmin = eps(T)
    backup_pose = copy(pose)
    
    #region MAINLOOP
    while driver_state.step < driver.max_steps

        # calculate current scaling factor and step size
        γ = min(γ, driver.max_displacement)
        driver_state.stepsize = γ / driver_state.max_force[1]
        # driver_state.stepsize = 1e-1

        # update coordinates
        for atom_index in 1:pose.state.size
            t = driver_state.stepsize .* pose.state.f[:, atom_index]
            pose.state.x[:, atom_index] = pose.state.x[:, atom_index] .+ t
        end

        # Sync pose since the "eval!" method might require internal coordinates
        ProtoSyn.request_c2i!(pose.state)
        sync!(pose)

        # Calculate new energy and forces (make sure to reset forces)
        ProtoSyn.reset_forces!(pose.state)
        energy = driver.eval!(pose, update_forces_overwrite = true)
        driver_state.max_force = ProtoSyn.atmax(pose.state, :f)

        # Verify convergence
        driver_state.converged = driver_state.max_force[1] < driver.force_tolerance
        driver_state.stalled = γ < γmin

        if driver_state.converged || driver_state.stalled
            break
        end

        if energy >= backup_pose.state.e[:Total]
            γ *= γdec
            ProtoSyn.recoverfrom!(pose, backup_pose)
            driver_state.max_force = ProtoSyn.atmax(pose.state, :f)
        else
            backup_pose = ProtoSyn.recoverfrom!(backup_pose, pose)
            γ *= γinc
        end

        # update current step and call calback functions (if any)
        driver_state.step += 1
        driver.callback !== nothing && driver.callback(pose, driver_state)

    end # while
    #endregion

    driver_state.completed = true
    
    if driver_state.stalled
        println("Steepest descent converged to maximum machine precision.")
    end

    if driver_state.converged
        println("Steepest descent converged to minimum force tolerance ($(driver_state.max_force[1]) < $(driver.force_tolerance))")
    end
    
    return pose
end

# * Show -----------------------------------------------------------------------

function Base.show(io::IO, sdd::SteepestDescent, level_code::Opt{LevelCode} = nothing)
    init_level_code = level_code === nothing ? LevelCode() : level_code
    init_lead       = ProtoSyn.get_lead(level_code)
    init_inner_lead = ProtoSyn.get_inner_lead(level_code)

    println(io, init_lead*"⚒  Steepest Descent Driver:")

    level_code = vcat(init_level_code, 3)
    lead       = ProtoSyn.get_lead(level_code)
    inner_lead = ProtoSyn.get_inner_lead(level_code)

    if typeof(sdd.eval!) === ProtoSyn.Calculators.EnergyFunction
        println(io, lead*" ●  Evaluator:")
        Base.show(io, sdd.eval!, vcat(level_code, 4))
    else
        println(io, lead*" ●  Evaluator: $(sdd.eval!)")
    end

    println(io, init_inner_lead*level_code.code_table[1])

    if sdd.callback !== nothing
        println(io, lead*" ● Callback:")
        Base.show(io, sdd.callback, vcat(level_code, 4))
    else
        println(io, lead*" ○  Callback: Not set")
    end

    level_code = vcat(init_level_code, 4)
    lead       = ProtoSyn.get_lead(level_code)
    inner_lead = ProtoSyn.get_inner_lead(level_code)

    println(io, init_inner_lead*level_code.code_table[1])

    println(io, lead*" ● Settings:")
    println(io, inner_lead*"       Max steps: $(sdd.max_steps)")
    println(io, inner_lead*" Force tolerance: $(sdd.force_tolerance)")
    println(io, inner_lead*"Max displacement: $(sdd.max_displacement)")
end
