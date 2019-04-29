module MD

# using ProtoSyn
using ..Common

mutable struct DriverOptions{F<:Function}
    evaluator!::F
    temperature::Float64
    n_steps::Int
    timestep::Float64
    tcoupling::Float64
end

function get_options(args)
    options = Dict{Symbol,Symbol}()
    for kwarg in args
        if Meta.isexpr(kwarg, :(=))
            key,val = kwarg.args
            options[key] = val
        else
            throw(ArgumentError("non-keyword argument like option '$kwarg'"))
        end
    end
    options
end





macro md(args...)
    
    options = get_options(args)

ex = quote
function(state::Common.State, driver::DriverOptions)
    
    n_atoms = state.size
    step::Int = 0       # current step
    Ti = 0.0            # instant temperature

    mass = 6.0 * ones(n_atoms)

    Δt = driver.timestep
    half_Δt = 0.5*Δt

    coords = state.xyz
    forces = state.forces
    #println(stderr, pointer_from_objref(coords))
    #println(stderr, pointer_from_objref(state.xyz))

    #velocs = zeros(n_atoms, 3)

    forces_t = copy(forces) # forces at time t

    evaluator! = driver.evaluator!
    
    energy = evaluator!(state, true)
    Ti = 0.0
    sum_sq = 0.0
    λ = 0.0
    Ndf::Int = 3*n_atoms - 6
    boltzmann = 0.0083145112119   # (kJ/(mol K))
    total_mass = sum(mass)

    # TEMPERATURE ASSIGNMENT
    #  1. initial temperatures
    λ = driver.temperature*boltzmann
    velocs = (*).(repeat(λ*mass,1,3), randn(n_atoms, 3))
    
    #  2. remove COM motion
    velocs .-= sum((*).(velocs, repeat(mass,1,3)), dims=1)
    
    #  3. scale to desired temperature
    Ti = 0.0
    @inbounds for i = 1:n_atoms
        sum_sq  = velocs[i,1]*velocs[i,1]
        sum_sq += velocs[i,2]*velocs[i,2]
        sum_sq += velocs[i,3]*velocs[i,3]
        Ti += sum_sq*mass[i]
    end
    Ti /= (Ndf*boltzmann)
    velocs .*= sqrt(driver.temperature/Ti)
    


    while step <= driver.n_steps
        
        # calculate instant temperature
        Ti = 0.0
        @inbounds for i = 1:n_atoms
            sum_sq  = velocs[i,1]*velocs[i,1]
            sum_sq += velocs[i,2]*velocs[i,2]
            sum_sq += velocs[i,3]*velocs[i,3]
            Ti += sum_sq*mass[i]
        end
        Ti /= (Ndf*boltzmann)
        
        # =============== START CONFIGURABLE SECTION ================
        #                 THERMOSTATS
        $(
        if get(options, :thermostat, :none) == :berendsen
            :(λ = sqrt(1.0 + Δt*(driver.temperature/Ti- 1.0)/driver.tcoupling))
        
        elseif get(options, :thermostat, :none) == :vrescale
            :(λ = sqrt(driver.temperature/Ti))
        end
        )
        
        $(if haskey(options, :thermostat)
            # :(velocs .*= max(0.85, min(1.25, λ)))
            :(velocs .*= λ)
        end)
        # ================ END CONFIGURABLE SECTION =================
        
        # println(stderr, step, " λ=", λ, " Ti=", Ti)
        
        # if step%50 == 0
        #     print(Print.as_xyz(state, metadata))
        # end

        # integrate equation of motion (Verlet velocity)
        @inbounds for i = 1:n_atoms
            f = half_Δt/mass[i]
            for k=1:3
                coords[i, k] += Δt * (f*forces[i,k] + velocs[i,k])
                forces_t[i,k] = forces[i,k]
            end
        end

        # calculate forces at t + Δt
        energy = evaluator!(state, true)
        
        # update velocities
        @inbounds for i = 1:n_atoms
            f = half_Δt/mass[i]
            @inbounds for k=1:3
                velocs[i,k] += f*(forces[i,k] + forces_t[i,k])
            end
        end

        # =============== START CONFIGURABLE SECTION ================
        # remove linear/angular momentum ----------------------------
        #if (this.momRemover && (step%20 == 0)) {
        #  this.momRemover.remove(coords, velocities, masses);
        #}

        # $(
        # if get(options, :com, :none) == :angular
        #     quote
        #         for i = 1:n_atoms
        #             m = masses[i]
        #             vi = @view velocs[i,:]
        #             xi = @view coords[i,:]
        #
        #             @. gx += m*xi               # group COM
        #             @. gp += m*vi               # group linear momentum
        #             @. gj += m*cross(xi, vi)    # group angular momentum
        #             gi .+= m*(xi * xi')         # group inertia tensor
        #
        #         end
        #
        #         gv[:] .= gp ./ total_mass # group velocity
        #         gx ./= total_mass # calculate group center of mass
        #
        #         # subtract COM contribution to the angular momentum
        #         gj[:] .-= total_mass .* cross(gx, gv)
        #
        #         # Subtract the center of mass contribution from the inertia tensor
        #         Icm[:] .= total_mass .* (gx * gx')
        #         gi .-= ICM
        #
        #         tmp[1,1] = gi[2,2]+gi[3,3]
        #         tmp[2,2] = gi[1,1]+gi[3,3]
        #         tmp[3,3] = gi[1,1]+gi[2,2]
        #         tmp[2,1] = tmp[1,2] = -gi[1,2]
        #         tmp[3,2] = tmp[2,3] = -gi[2,3]
        #         tmp[3,1] = tmp[1,3] = -gi[1,3]
        #         rfac = (tmp[1,1]+tmp[2,2]+tmp[3,3])/3
        #         tmp .*= rfac
        #
        #     end
        # end
        # ) # end interpolation section
        # ================ END CONFIGURABLE SECTION =================

        
        step += 1
    end

end # end fcn declaration
end # end quote
    ex
end


run! = @md thermostat=vrescale


end