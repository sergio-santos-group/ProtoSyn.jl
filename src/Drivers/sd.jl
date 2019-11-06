
#------------------------------------------------
Base.@kwdef mutable struct SteepestDescentState <: IDriverState
    step::Int64        = 0
    converged::Bool    = false
    completed::Bool    = false
    stalled::Bool      = false

    stepsize::Float64 =  1.0
    max_force::Float64 = -1.0
end


#------------------------------------------------
Base.@kwdef mutable struct SteepestDescent{F <: Function} <: IDriver
    eval!::F
    max_steps::Int = 0
    pairlist_freq::Int = 0

    force_tolerance::Float64 = 1.0
    max_displacement::Float64 = 0.1
end


#------------------------------------------------
function (driver::SteepestDescent)(cb::Opt{F}, state::State) where {F <: Function}

    if state.forces === nothing 
        state.forces = zeros(state.size, 3)
    end

    if state.energy === nothing 
        state.energy = ProtoSyn.Energy()
    end


    # Evaluate initial energy and forces
    #   start by calculating nonbonded lists:
    #   assuming that each atomic displacement is limited
    #   to 'max_displacement' at each step, then any
    #   2 atoms, if moving in oposite directions, can be
    #   displaced 2*max_displacement from each other.
    #   Therefore, the maximum total displacement during
    #   'pairlist_freq' steps is 2*max_displacement*pairlist_freq.
    #   Hence, the nblist cutoff is (cutoff + 2*max_displacement*pairlist_freq).
    
    if state.pairlist !== nothing
        state.pairlist.buffer  = 2.0 * driver.max_displacement * driver.pairlist_freq
        nbupdate!(state)
    end
    
    energy = driver.eval!(state, true)

    # instantiate a new DriverState object.
    # By default, no optimization step has yet been taken
    # apart from calculating the energy and forces for the
    # input state
    driver_state = SteepestDescentState()
    driver_state.max_force = max(state, :forces)
    
    # if max force is already below the requested force tolerance,
    # simply return the current driver state
    driver_state.converged = driver_state.max_force < driver.force_tolerance

    if driver_state.converged
        return driver_state
    end
    
    # Initial callbacks
    cb !== nothing && cb(state, driver_state)
    

    # Initialize γ and backup copy
    γ = 1.0
    backup_state = State(state)
    
    #region MAINLOOP
    while driver_state.step < driver.max_steps
        
        # save current state to backup
        copy!(backup_state, state)

        # calculate current scaling factor and step size
        γ = min(γ, driver.max_displacement)
        driver_state.stepsize = γ / driver_state.max_force

        # update coordinates
        @. state.coords += driver_state.stepsize * state.forces

        # update nonbonded lists if required
        # TODO: check nblists 
        if (driver.pairlist_freq > 0) &&
           (driver_state.step > 0) &&
           (driver_state.step % driver.pairlist_freq == 0)
           nbupdate!(state)
        end

        # Calculate new energy and forces
        # (make sure to reset forces)
        fill!(state.forces, 0.0)
        energy = driver.eval!(state, true)
        driver_state.max_force = max(state, :forces)

        # Verify convergence
        driver_state.converged = driver_state.max_force < driver.force_tolerance
        driver_state.stalled = γ < eps()

        if driver_state.converged || driver_state.stalled
            break
        end
        
        # Update gamma
        if energy >= backup_state.energy.total
            γ *= 0.50
            copy!(state, backup_state)
            driver_state.max_force = max(state, :forces)
        else
            γ *= 1.05
        end

        # update current step and call calback functions (if any)
        driver_state.step += 1
        cb !== nothing && cb(state, driver_state)

    end # while
    #endregion

    driver_state.completed = true
    
    driver_state
end

(driver::SteepestDescent)(s::State) = driver(nothing, s)


