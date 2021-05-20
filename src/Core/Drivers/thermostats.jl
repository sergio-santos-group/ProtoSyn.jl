# abstract type IThermostat end


# struct BerendsenThermostat <: IThermostat
#     temperature::Float64
#     coupling::Float64
#     timestep::Float64
# end

# @inline function (th::BerendsenThermostat)(state::State, t::Float64)
#     λ = sqrt(1.0 + th.timestep*(th.temperature/t - 1.0) / th.coupling)
#     @. state.velocs *= λ
#     state
# end


# struct VRescaleThermostat <: IThermostat
#     temperature::Float64
# end

# @inline function (th::VRescaleThermostat)(state::State, t::Float64)
#     λ = sqrt(th.temperature/t)
#     @. state.velocs *= λ
#     state
# end

"""
    # TODO
"""
function get_constant_temperature(init_temp::T) where {T <: AbstractFloat}
    return function constant_temperature(step::Int)
        return init_temp
    end
end

"""
    # TODO
"""
function get_linear_quench(init_temp::T, n_steps::Int) where {T <: AbstractFloat}
    return function linear_quench(step::Int)
        return -((step * init_temp)/n_steps) + init_temp
    end
end

"""
    # TODO
"""
function get_sigmoid_quench(init_temp::T, quench_rate::T, min_temp::T) where {T <: AbstractFloat}
    return function sigmoid_quench(step::Int)
        return ((1/(1+(exp(step/(quench_rate+(quench_rate/10))))))*2*init_temp-init_temp)+init_temp
    end
end