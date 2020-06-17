# using ..Calculators

# #------------------------------------------------
# Base.@kwdef mutable struct SteepestDescentState{T<:AbstractFloat} <: IDriverState
# # Base.@kwdef mutable struct SteepestDescentState <: IDriverState
#     step::Int          = 0
#     converged::Bool    = false
#     completed::Bool    = false
#     stalled::Bool      = false

#     stepsize::T =  T(1)
#     max_force::T = T(-1)
#     # stepsize::Float64 =  1.0
#     # max_force::Float64 = -1.0
# end


# #------------------------------------------------
# Base.@kwdef mutable struct SteepestDescent{F <: Function} <: IDriver
#     eval!::F
#     max_steps::Int = 0
#     pairlist_freq::Int = 0

#     force_tolerance::Float64 = 1.0
#     max_displacement::Float64 = 0.1
# end


# #------------------------------------------------
# function (driver::SteepestDescent)(cb::Opt{F}, state::Calculators.State) where {F <: Function}

#     # Evaluate initial energy and forces
#     #   start by calculating nonbonded lists:
#     #   assuming that each atomic displacement is limited
#     #   to 'max_displacement' at each step, then any
#     #   2 atoms, if moving in oposite directions, can be
#     #   displaced 2*max_displacement from each other.
#     #   Therefore, the maximum total displacement during
#     #   'pairlist_freq' steps is 2*max_displacement*pairlist_freq.
#     #   Hence, the nblist cutoff is (cutoff + 2*max_displacement*pairlist_freq).
    
#     #if state.pairlist !== nothing
#     #    state.pairlist.buffer  = 2.0 * driver.max_displacement * driver.pairlist_freq
#     #    nbupdate!(state)
#     #end
    
#     energy = driver.eval!(state, true)

#     # instantiate a new DriverState object.
#     # By default, no optimization step has yet been taken
#     # apart from calculating the energy and forces for the
#     # input state
#     driver_state = SteepestDescentState()
#     driver_state.max_force = atmax(state, :f)
    
#     # if max force is already below the requested force tolerance,
#     # simply return the current driver state
#     driver_state.converged = driver_state.max_force < driver.force_tolerance

#     if driver_state.converged
#         return driver_state
#     end
    
#     # Initial callbacks
#     cb !== nothing && cb(state, driver_state)
    

#     # Initialize γ and backup copy
#     γ = 1.0
#     backup_state = Calculators.State(state)
    
#     #region MAINLOOP
#     while driver_state.step < driver.max_steps
        
#         # save current state to backup
#         copy!(backup_state, state)

#         # calculate current scaling factor and step size
#         γ = min(γ, driver.max_displacement)
#         driver_state.stepsize = γ / driver_state.max_force

#         # update coordinates
#         # @. state.coords += driver_state.stepsize * state.forces
#         @. state.x += driver_state.stepsize * state.f

#         # update nonbonded lists if required
#         # TODO: check nblists 
#         # if (driver.pairlist_freq > 0) &&
#         #    (driver_state.step > 0) &&
#         #    (driver_state.step % driver.pairlist_freq == 0)
#         #    nbupdate!(state)
#         # end

#         # Calculate new energy and forces
#         # (make sure to reset forces)
#         #fill!(state.forces, 0.0)
#         energy = driver.eval!(state, true)
#         driver_state.max_force = atmax(state, :f)

#         # Verify convergence
#         driver_state.converged = driver_state.max_force < driver.force_tolerance
#         driver_state.stalled = γ < eps()

#         if driver_state.converged || driver_state.stalled
#             break
#         end
        
#         # Update gamma
#         # if energy >= backup_state.energy.total
#         if energy >= Calculators.energy(backup_state)
#             γ *= 0.50
#             copy!(state, backup_state)
#             driver_state.max_force = atmax(state, :f)
#         else
#             γ *= 1.05
#         end

#         # update current step and call calback functions (if any)
#         driver_state.step += 1
#         cb !== nothing && cb(state, driver_state)

#     end # while
#     #endregion

#     driver_state.completed = true
    
#     driver_state
# end

# (driver::SteepestDescent)(s::Calculators.State) = driver(nothing, s)


# function atmax(state::Calculators.State{T}, comp::Symbol) where T
#     m = T(0)
#     x = getproperty(state, comp)
#     for i = 1:state.size
#         f = x[i,1]^2 + x[i,2]^2 + x[i,3]^2
#         if f > m
#             m = f
#         end
#     end
#     return sqrt(m)
# end




using ..Calculators

#------------------------------------------------
Base.@kwdef mutable struct SteepestDescentState{T<:AbstractFloat} <: IDriverState
    step::Int          = 0
    converged::Bool    = false
    completed::Bool    = false
    stalled::Bool      = false

    stepsize::T =  T(1)
    max_force::Tuple{T,Int} = (T(-1),0)
end


#------------------------------------------------
Base.@kwdef mutable struct SteepestDescent{T<:AbstractFloat} <: IDriver
    eval!::Function
    max_steps::Int = 0
    pairlist_freq::Int = 0

    force_tolerance::T = T(1.0)
    max_displacement::T = T(0.1)
end


#------------------------------------------------
function (driver::SteepestDescent{T})(cb::Opt{F}, state::Calculators.State{T}) where {T, F <: Function}

    # Evaluate initial energy and forces
    #   start by calculating nonbonded lists:
    #   assuming that each atomic displacement is limited
    #   to 'max_displacement' at each step, then any
    #   2 atoms, if moving in oposite directions, can be
    #   displaced 2*max_displacement from each other.
    #   Therefore, the maximum total displacement during
    #   'pairlist_freq' steps is 2*max_displacement*pairlist_freq.
    #   Hence, the nblist cutoff is (cutoff + 2*max_displacement*pairlist_freq).
    
    #if state.pairlist !== nothing
    #    state.pairlist.buffer  = 2.0 * driver.max_displacement * driver.pairlist_freq
    #    nbupdate!(state)
    #end
    
    energy = driver.eval!(state, true)

    # instantiate a new DriverState object.
    # By default, no optimization step has yet been taken
    # apart from calculating the energy and forces for the
    # input state
    driver_state = SteepestDescentState{T}()
    driver_state.max_force = atmax(state, :f)
    
    # if max force is already below the requested force tolerance,
    # simply return the current driver state
    driver_state.converged = driver_state.max_force[1] < driver.force_tolerance

    if driver_state.converged
        return driver_state
    end
    
    # Initial callbacks
    cb !== nothing && cb(state, driver_state)
    

    # Initialize γ and backup copy
    γ = T(1.0)
    γmin = eps(T)
    γinc = T(1.05)
    γdec = T(0.50)
    backup_state = Calculators.State(state)
    
    #region MAINLOOP
    while driver_state.step < driver.max_steps
        
        # save current state to backup
        copy!(backup_state, state)

        # calculate current scaling factor and step size
        γ = min(γ, driver.max_displacement)
        driver_state.stepsize = γ / driver_state.max_force[1]

        # update coordinates
        # @. state.coords += driver_state.stepsize * state.forces
        @. state.x += driver_state.stepsize * state.f

        # update nonbonded lists if required
        # TODO: check nblists 
        # if (driver.pairlist_freq > 0) &&
        #    (driver_state.step > 0) &&
        #    (driver_state.step % driver.pairlist_freq == 0)
        #    nbupdate!(state)
        # end

        # Calculate new energy and forces
        # (make sure to reset forces)
        energy = driver.eval!(state, true)
        driver_state.max_force = atmax(state, :f)

        # Verify convergence
        driver_state.converged = driver_state.max_force[1] < driver.force_tolerance
        driver_state.stalled = γ < γmin

        if driver_state.converged || driver_state.stalled
            break
        end
        
        # Update gamma
        if energy >= Calculators.energy(backup_state)
            γ *= γdec
            state, backup_state = backup_state, state
            # copy!(state, backup_state)
            driver_state.max_force = atmax(state, :f)
        else
            γ *= γinc
        end

        # update current step and call calback functions (if any)
        driver_state.step += 1
        cb !== nothing && cb(state, driver_state)

    end # while
    #endregion

    driver_state.completed = true
    
    driver_state
end

(driver::SteepestDescent)(s::Calculators.State) = driver(nothing, s)

