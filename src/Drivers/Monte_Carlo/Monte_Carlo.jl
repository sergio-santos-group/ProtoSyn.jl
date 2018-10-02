module Monte_Carlo

using ..Common
using Printf

struct ConfigParameters

    nsteps::Int64
    temperature::Float64

end

#TO DO:
# 1) DOCUMENT THE FUNCTION
function load(p::Dict{String, Any})::ConfigParameters

    return ConfigParameters(p["nsteps"], p["temperature"])

end

#TO DO:
# 1) DOCUMENT THE FUNCTION
function run!(
    state::Common.State,
    sampler!::Function,
    evaluator!::Function,
    params::ConfigParameters;
    ostream::IO = stdout,
    callback::Union{Function, Nothing} = nothing)

    xyz0 = copy(state.xyz)
    ene0 = evaluator!(state, false)

    step::Int64 = 0
    while step < params.nsteps
        step += 1

        #Generate new state and evaluate its energy
        sampler!(state)
        ene1 = evaluator!(state, false)

        if ene1 < ene0 || rand() < exp(-(ene1 - ene0) / params.temperature)
            ene0 = ene1
            xyz0[:] = state.xyz
            if callback != nothing
                callback(state, step)
            end
            write(ostream, @sprintf "(MC) Step: %4d | Energy: %9.4f\n" step ene1)
        else
            state.xyz[:] = xyz0
        end
    end
end

end