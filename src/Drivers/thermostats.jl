abstract type IThermostat end


struct BerendsenThermostat <: IThermostat
    temperature::Float64
    coupling::Float64
    timestep::Float64
end

@inline function (th::BerendsenThermostat)(state::State, t::Float64)
    λ = sqrt(1.0 + th.timestep*(th.temperature/t - 1.0) / th.coupling)
    @. state.velocs *= λ
    state
end


struct VRescaleThermostat <: IThermostat
    temperature::Float64
end

@inline function (th::VRescaleThermostat)(state::State, t::Float64)
    λ = sqrt(th.temperature/t)
    @. state.velocs *= λ
    state
end