using ProtoSyn
using ProtoSyn.Calculators: EnergyFunction

Base.@kwdef mutable struct SteepestDescentState{T <: AbstractFloat} <: DriverState
    step::Int          = 0
    converged::Bool    = false
    completed::Bool    = false
    stalled::Bool      = false

    stepsize::T        =  T(1)
    max_force::Tuple{T,Int} = (T(-1),0)
end


mutable struct SteepestDescent{T <: AbstractFloat} <: Driver
    eval!::Union{Function, EnergyFunction}
    callback::Opt{Callback}
    max_steps::Int
    force_tolerance::T
    max_displacement::T
end


SteepestDescent(::Type{T}) where {T <: AbstractFloat} = begin
    zero_eval = (pose::Pose; update_forces = False) -> return 0.0
    return SteepestDescent{T}(zero_eval, 0, T(1), T(0.1))
end
SteepestDescent() = begin
    return SteepestDescent(ProtoSyn.Units.defaultFloat)
end

function Base.show(io::IO, sd::SteepestDescent)
    println(io, "Steepest Descent Driver:")
    println(io, " $(sd.eval!)")
    println(io, " Max Steps: $(sd.max_steps)")
    println(io, " Force Tolerance: $(sd.force_tolerance)")
    println(io, " Max displacement: $(sd.max_displacement)")
end


function (driver::SteepestDescent{T})(pose::Pose) where {T}
  
    energy = driver.eval!(pose, update_forces = true)

    # Onstantiate a new DriverState object.
    # By default, no optimization step has yet been taken
    # apart from calculating the energy and forces for the
    # input state
    driver_state = SteepestDescentState{T}()
    driver_state.max_force = ProtoSyn.atmax(pose.state, :f)
    
    # If max force is already below the requested force tolerance,
    # simply return the current driver state
    driver_state.converged = driver_state.max_force[1] < driver.force_tolerance

    if driver_state.converged
        return driver_state
    end
    
    # Initial callbacks
    driver.callback !== nothing && driver.callback(pose, driver_state)
    

    # Initialize γ and backup copy
    γ    = T(1.0)
    γmin = eps(T)
    γinc = T(1.05)
    γdec = T(0.50)
    # backup_pose = copy(pose)
    
    #region MAINLOOP
    while driver_state.step < driver.max_steps
        
        # save current state to backup
        backup_pose = copy(pose)

        # calculate current scaling factor and step size
        γ = min(γ, driver.max_displacement)
        # driver_state.stepsize = γ / driver_state.max_force[1]
        driver_state.stepsize = 1e-1

        # update coordinates
        for atom_index in 1:pose.state.size
            t = driver_state.stepsize .* pose.state.f[:, atom_index]
            pose.state.x[:, atom_index] = pose.state.x[:, atom_index] .- t
        end

        # Calculate new energy and forces
        # (make sure to reset forces)
        ProtoSyn.reset_forces!(pose.state)
        energy = driver.eval!(pose, update_forces = true)
        driver_state.max_force = ProtoSyn.atmax(pose.state, :f)

        # Verify convergence
        driver_state.converged = driver_state.max_force[1] < driver.force_tolerance
        driver_state.stalled = γ < γmin

        if driver_state.converged || driver_state.stalled
            break
        end
        
        # Update gamma (BUG HERE !!!!!!!!!!)
        # if energy >= driver.eval!(backup_pose)
        #     γ *= γdec
        #     pose, backup_pose = backup_pose, pose
        #     driver_state.max_force = ProtoSyn.atmax(pose.state, :f)
        # else
        #     γ *= γinc
        # end

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
    driver_state
end