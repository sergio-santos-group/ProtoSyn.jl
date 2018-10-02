module Steepest_Descent

using ..Common
using Printf

struct SDConfigParameters
    
    nsteps::Int64
    log_freq::Int64
    fTol::Float64
    max_step::Float64

end


#TO DO:
# 1) DOCUMENT THE FUNCTION
function load(p::Dict{String, Any})::SDConfigParameters

    return SDConfigParameters(p["nsteps"], p["log_freq"], p["fTol"], p["max_step"])

end

#TO DO:
# 1) GAMMA SHOULD BE UPDATED BASED ON ENERGY (?)
# 2) DOCUMENT THE FUNCTION
function run!(
    state::Common.State,
    evaluator!::Function,
    params::SDConfigParameters;
    ostream::IO = stdout,
    callback::Union{Function, Nothing} = nothing)

    function get_max_force(f::Array{Float64, 2})
        return sqrt(maximum(sum(f.^2, dims = 2)))
    end

    # Get total number of atoms
    natoms = size(state.xyz, 1)

    # Evaluate initial energy and forces
    energy_old::Float64 = evaluator!(state, true)
    max_force = get_max_force(state.forces)
    old_max_force::Float64 = max_force
    
    
    gamma::Float64 = 1e-5
    step::Int64 = 0
    if callback != nothing
        callback(state, step)
    end
    while step < params.nsteps
        step += 1

        #Update system coordinates
        gamma = min(gamma, params.max_step)
        dx = (gamma / get_max_force(state.forces)) * state.forces
        state.xyz += dx

        #Calculate new energy and forces
        state.forces = zeros(size(state.xyz, 1), 3)
        energy = evaluator!(state, true)
        max_force = get_max_force(state.forces)

        #Callback, if present
        if step % params.log_freq == 0
            if callback != nothing
                callback(state, step)
            end
            write(ostream, @sprintf "(SD) Step: %4d | Energy: %9.4f | maxForce: %9.4f\n" step energy max_force)
        end

        #If max force variance is bellow fTol threshold, convergence was achieved: exit function
        if abs(old_max_force - max_force) < params.fTol
            if callback != nothing
                callback(state, step)
            end
            write(ostream, "Achieved convergence (below fTol $(params.fTol)) in $step steps.\n")
            break
        end

        #Update gamma constant
        if max_force > old_max_force
            gamma *= 0.50
        else
            gamma *= 1.05
        end

        #Update current force to new cycle
        old_max_force = max_force
    end
end

end