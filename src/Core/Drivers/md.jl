using Base.Cartesian
using ..ProtoSyn: @cross
using LinearAlgebra: mul!



#------------------------------------------------
Base.@kwdef mutable struct MolecularDynamicsState <: IDriverState
    step::Int64        = 0
    converged::Bool    = false
    completed::Bool    = false
    stalled::Bool      = false

    temperature::Float64 =  -1.0    # instant temperature
end



#------------------------------------------------
Base.@kwdef mutable struct MolecularDynamics{F <: Function, T<:IThermostat} <: IDriver
    eval!::F
    max_steps::Int = 0
    pairlist_freq::Int = 0

    #temperature::Float64
    #tcoupling::Float64
    masses::Vector{Float64}
    timestep::Float64
    thermostat!::T

    remove_com_freq::Int
    remove_com_mode::Symbol

end



@inline function temperature(state::State, masses::Vector{Float64})
    t = 0.0
    v = state.velocs
    @inbounds for i = 1:state.size
        t += masses[i] * (v[i,1]^2 + v[i,2]^2 + v[i,3]^2)
    end
    t / ((3*state.size - 6) * ProtoSyn.Boltzmann)
end



function (driver::MolecularDynamics)(cb::Opt{F}, state::State) where {F <: Function}
    if state.forces === nothing 
        state.forces = zeros(state.size, 3)
    end
    if state.energy === nothing 
        state.energy = ProtoSyn.Energy()
    end
    
    T = driver.thermostat!.temperature
    # mean speed (m/s) for a free Hydrogen atom
    #   R is in kJ mol-1 K-1
    # divide by 1e3 to convert to nm/ps
    #  mean_speed = sqrt(8*ProtoSyn.Boltzmann*T*1e6/(π*m)) / 1000
    m = 1.0
    mean_speed = sqrt(8*ProtoSyn.Boltzmann*T/(π*m))
    
    # the maximum displacement of a particle travelling
    # at mean_speed during "pairlist_freq" steps is
    #  max_disp = mean_speed * pairlist_freq * timestep.
    # Hence, any two atoms travelling in opposite directions
    # can become apart by 2*max_disp nm. The verlet list
    # buffer is thus 2*max_disp.
    if state.pairlist !== nothing
        max_disp = mean_speed * driver.pairlist_freq * driver.timestep
        println("mean_speed = $mean_speed  buffer = $(2*max_disp)")
        state.pairlist.buffer = 2.0 * max_disp
        nbupdate!(state)
    end
    
    dt = driver.timestep
    half_dt_over_mass = (0.5*dt) ./ driver.masses
    
    driver_state = MolecularDynamicsState()
    
    energy = driver.eval!(state, true)
    
    # check if any velocities already exist. If not,
    # generate velocities for the requested temperature
    if state.velocs === nothing
        
        # determine target temperature and generate random
        # velocities from a boltzmann distribution at T0
        T0 = driver.thermostat!.temperature 
        λ = T0 * ProtoSyn.Boltzmann
        state.velocs = repeat(λ*driver.masses,1,3) .* randn(state.size, 3)

        # remove center of mass motion
        removecom!(state, driver.masses, :linear)

        #  and rescale to desired temperature
        driver_state.temperature = temperature(state, driver.masses)
        state.velocs .*= sqrt(T0/driver_state.temperature)
        
    else
        # otherwise, simply recalculate the current
        # instantaneous temperature
        driver_state.temperature = temperature(state, driver.masses)
    end

    
    velocs = state.velocs   # velocities at time t+dt
    coords = state.coords   # coordinates at time t+dt
    forces = state.forces   # forces at time t+dt
    forces_t = copy(forces) # forces at time t
    
    # call callback function (if provided)
    cb !== nothing && cb(state, driver_state)


    # MAIN INTEGRATION LOOP
    while driver_state.step < driver.max_steps
        
        driver_state.step += 1

        # determine instantaneous temperature
        driver_state.temperature = temperature(state, driver.masses)

        # rescale velocities using the provided thermostat
        driver.thermostat!(state, driver_state.temperature)

        # integrate equation of motion (Verlet velocity)
        @inbounds for i = 1:state.size
            f = half_dt_over_mass[i]
            @nexprs 3 u -> coords[i, u] += dt * (f*forces[i, u] + velocs[i, u])
            @nexprs 3 u -> forces_t[i, u] = forces[i, u]
        end

        # calculate forces at t + dt
        energy = driver.eval!(state, true)

        # update velocities at t + dt
        @inbounds for i = 1:state.size
            f = half_dt_over_mass[i]
            @nexprs 3 u -> velocs[i, u] += f*(forces[i, u] + forces_t[i, u])
        end

        # remove COM motion
        if (driver.remove_com_freq > 0) &&
            (driver_state.step % driver.remove_com_freq == 0)
            removecom!(state, driver.masses, driver.remove_com_mode)
        end
        
        # update nonbonded lists if required
        if (driver.pairlist_freq > 0) &&
            (driver_state.step > 0) &&
            (driver_state.step % driver.pairlist_freq == 0)
            nbupdate!(state)
         end

        # call callback function (if provided)
        cb !== nothing && cb(state, driver_state)

    end # while MAIN INTEGRATION LOOP

    driver_state.completed = true
    driver_state

end
(driver::MolecularDynamics)(s::State) = driver(nothing, s)



function removecom!(state::State, masses::Vector{Float64}, mode::Symbol)
    
    velocs = state.velocs
    coords = state.coords
    
    # linear momentum remover --------------------------------------------------
    if mode==:linear
        #println("removecom linear")
        mtot = 0.0
        @nexprs 3 u -> vcom_u = 0.0
        @inbounds for i=1:state.size
            mass = masses[i]
            mtot += mass
            @nexprs 3 u -> vcom_u += mass * velocs[i,u]
        end

        @nexprs 3 u -> vcom_u /= mtot
        @inbounds for i=1:state.size
            @nexprs 3 u -> velocs[i,u] -= vcom_u
        end

    # angular momentum remover -------------------------------------------------
    elseif mode == :angular
        #println("removecom angular")
        
        tmp = zeros(3,3)
        gw  = zeros(3)
        gj  = zeros(3)

        @nexprs 3 u -> gx_u = 0.0
        @nexprs 3 u -> gp_u = 0.0
        @nexprs 3 u -> gj_u = 0.0
        @nexprs 3 k -> @nexprs 3 u -> gi_u_k = 0.0
        
        mtot = 0.0
        @inbounds for i = 1:state.size
            mass = masses[i]
            mtot += mass
            # center of mass
            @nexprs 3 u -> gx_u += mass * coords[i, u]
            
            # linear momentum
            @nexprs 3 u -> gp_u += mass * velocs[i, u]
            
            # angular momentum
            @cross u xv_u coords[i,u] velocs[i,u]
            @nexprs 3 u -> gj_u += mass*xv_u
            
            # group inertia tensor
            @nexprs 3 k -> @nexprs 3 u -> gi_u_k += mass * coords[i,u] * coords[i,k]
        end

        # group velocity
        @nexprs 3 u -> gv_u = gp_u / mtot
        
        # group center of mass
        @nexprs 3 u -> gx_u /= mtot

        # subtract COM contribution to the angular momentum
        @cross u xv_u gx_u gv_u
        @nexprs 3 u -> gj[u] = gj_u - mtot * xv_u

        @nexprs 3 k -> @nexprs 3 u -> Icm_u_k = mtot * gx_u * gx_k
        @nexprs 3 k -> @nexprs 3 u ->  gi_u_k -= Icm_u_k

        tmp[1,1] = gi_2_2 + gi_3_3
        tmp[2,2] = gi_1_1 + gi_3_3
        tmp[3,3] = gi_1_1 + gi_2_2
        tmp[2,1] = tmp[1,2] = -gi_1_2
        tmp[3,2] = tmp[2,3] = -gi_2_3
        tmp[3,1] = tmp[1,3] = -gi_1_3
        fac = 3.0 / (tmp[1,1] + tmp[2,2] + tmp[3,3])
        tmp .*= fac
        itmp = inv(tmp)
        itmp .*= fac
        mul!(gw, itmp, gj)

        @inbounds for i=1:state.size
            @nexprs 3 u -> velocs[i,u] -= gv_u
            @nexprs 3 u -> dx_u = coords[i,u] - gx_u
            @cross u dv_u gw[u] dx_u
            @nexprs 3 u -> velocs[i,u] -= dv_u
        end
    end

    state
end




# macro md(args...)

#     options = get_options(args)

#     ex = quote
#     function(state::Common.State, driver::DriverOptions)

#         n_atoms = state.size
#         step::Int = 0       # current step
#         Ti = 0.0            # instant temperature

#         mass = 6.0 * ones(n_atoms)

#         dt = driver.timestep
#         half_dt = 0.5*dt

#         coords = state.xyz
#         forces = state.forces
#         #println(stderr, pointer_from_objref(coords))
#         #println(stderr, pointer_from_objref(state.xyz))

#         #velocs = zeros(n_atoms, 3)

#         forces_t = copy(forces) # forces at time t

#         evaluator! = driver.evaluator!

#         energy = evaluator!(state, true)
#         Ti = 0.0
#         sum_sq = 0.0
#         λ = 0.0
#         Ndf::Int = 3*n_atoms - 6
#         boltzmann = 0.0083145112119   # (kJ/(mol K))
#         mtot = sum(mass)

#         # TEMPERATURE ASSIGNMENT
#         #  1. initial temperatures
#         λ = driver.temperature*boltzmann
#         velocs = (*).(repeat(λ*mass,1,3), randn(n_atoms, 3))

#         #  2. remove COM motion
#         velocs .-= sum((*).(velocs, repeat(mass,1,3)), dims=1)

#         #  3. scale to desired temperature
#         Ti = 0.0
#         @inbounds for i = 1:n_atoms
#             sum_sq  = velocs[i,1]*velocs[i,1]
#             sum_sq += velocs[i,2]*velocs[i,2]
#             sum_sq += velocs[i,3]*velocs[i,3]
#             Ti += sum_sq*mass[i]
#         end
#         Ti /= (Ndf*boltzmann)
#         velocs .*= sqrt(driver.temperature/Ti)



#         while step <= driver.n_steps

#             # calculate instant temperature
#             Ti = 0.0
#             @inbounds for i = 1:n_atoms
#                 sum_sq  = velocs[i,1]*velocs[i,1]
#                 sum_sq += velocs[i,2]*velocs[i,2]
#                 sum_sq += velocs[i,3]*velocs[i,3]
#                 Ti += sum_sq*mass[i]
#             end
#             Ti /= (Ndf*boltzmann)

#             # =============== START CONFIGURABLE SECTION ================
#             #                 THERMOSTATS
#             $(
#             if get(options, :thermostat, :none) == :berendsen
#                 :(λ = sqrt(1.0 + dt*(driver.temperature/Ti- 1.0)/driver.tcoupling))

#             elseif get(options, :thermostat, :none) == :vrescale
#                 :(λ = sqrt(driver.temperature/Ti))
#             end
#             )

#             $(if haskey(options, :thermostat)
#                 # :(velocs .*= max(0.85, min(1.25, λ)))
#                 :(velocs .*= λ)
#             end)
#             # ================ END CONFIGURABLE SECTION =================

#             # println(stderr, step, " λ=", λ, " Ti=", Ti)

#             # if step%50 == 0
#             #     print(Print.as_xyz(state, metadata))
#             # end

#             # integrate equation of motion (Verlet velocity)
#             @inbounds for i = 1:n_atoms
#                 f = half_dt/mass[i]
#                 for k=1:3
#                     coords[i, k] += dt * (f*forces[i,k] + velocs[i,k])
#                     forces_t[i,k] = forces[i,k]
#                 end
#                 Ti /= (Ndf*boltzmann)
#                 velocs .*= sqrt(driver.temperature/Ti)



#                 while step <= driver.n_steps

#                     # calculate instant temperature
#                     Ti = 0.0
#                     @inbounds for i = 1:n_atoms
#                         sum_sq  = velocs[i,1]*velocs[i,1]
#                         sum_sq += velocs[i,2]*velocs[i,2]
#                         sum_sq += velocs[i,3]*velocs[i,3]
#                         Ti += sum_sq*mass[i]
#                     end
#                     Ti /= (Ndf*boltzmann)

#                     # =============== START CONFIGURABLE SECTION ================
#                     #                 THERMOSTATS
#                     $(
#                     if get(options, :thermostat, :none) == :berendsen
#                         :(λ = sqrt(1.0 + dt*(driver.temperature/Ti- 1.0)/driver.tcoupling))

#                     elseif get(options, :thermostat, :none) == :vrescale
#                         :(λ = sqrt(driver.temperature/Ti))
#                     end
#                     )

#                     $(if haskey(options, :thermostat)
#                         # :(velocs .*= max(0.85, min(1.25, λ)))
#                         :(velocs .*= λ)
#                     end)
#                     # ================ END CONFIGURABLE SECTION =================

#                     # println(stderr, step, " λ=", λ, " Ti=", Ti)

#                     # if step%50 == 0
#                     #     print(Print.as_xyz(state, metadata))
#                     # end

#                     # integrate equation of motion (Verlet velocity)
#                     @inbounds for i = 1:n_atoms
#                         f = half_dt/mass[i]
#                         for k=1:3
#                             coords[i, k] += dt * (f*forces[i,k] + velocs[i,k])
#                             forces_t[i,k] = forces[i,k]
#                         end
#                     end

#                     # calculate forces at t + dt
#                     energy = evaluator!(state, true)

#                     # update velocities
#                     @inbounds for i = 1:n_atoms
#                         f = half_dt/mass[i]
#                         @inbounds for k=1:3
#                             velocs[i,k] += f*(forces[i,k] + forces_t[i,k])
#                         end
#                     end

#                     # =============== START CONFIGURABLE SECTION ================
#                     # remove linear/angular momentum ----------------------------
#                     #if (this.momRemover && (step%20 == 0)) {
#                     #  this.momRemover.remove(coords, velocities, masses);
#                     #}

#                     # $(
#                     # if get(options, :com, :none) == :angular
#                     #     quote
#                     #         for i = 1:n_atoms
#                     #             m = masses[i]
#                     #             vi = @view velocs[i,:]
#                     #             xi = @view coords[i,:]
#                     #
#                     #             @. gx += m*xi               # group COM
#                     #             @. gp += m*vi               # group linear momentum
#                     #             @. gj += m*cross(xi, vi)    # group angular momentum
#                     #             gi .+= m*(xi * xi')         # group inertia tensor
#                     #
#                     #         end
#                     #
#                     #         gv[:] .= gp ./ mtot # group velocity
#                     #         gx ./= mtot # calculate group center of mass
#                     #
#                     #         # subtract COM contribution to the angular momentum
#                     #         gj[:] .-= mtot .* cross(gx, gv)
#                     #
#                     #         # Subtract the center of mass contribution from the inertia tensor
#                     #         Icm[:] .= mtot .* (gx * gx')
#                     #         gi .-= ICM
#                     #
#                     #         tmp[1,1] = gi[2,2]+gi[3,3]
#                     #         tmp[2,2] = gi[1,1]+gi[3,3]
#                     #         tmp[3,3] = gi[1,1]+gi[2,2]
#                     #         tmp[2,1] = tmp[1,2] = -gi[1,2]
#                     #         tmp[3,2] = tmp[2,3] = -gi[2,3]
#                     #         tmp[3,1] = tmp[1,3] = -gi[1,3]
#                     #         rfac = (tmp[1,1]+tmp[2,2]+tmp[3,3])/3
#                     #         tmp .*= rfac
#                     #
#                     #     end
#                     # end
#                     # ) # end interpolation section
#                     # ================ END CONFIGURABLE SECTION =================


#                     step += 1
#                 end

#             end # end fcn declaration
#         end # end quote
#         ex
#     end
# end


# # run! = @md thermostat=vrescale


# end