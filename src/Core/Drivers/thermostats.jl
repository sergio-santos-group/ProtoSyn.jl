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
    get_constant_temperature(init_temp::T) where {T <: AbstractFloat}

Return a constant temperature thermostat function, set at temperature
`init_temp`.

# See also
[`get_linear_quench`](@ref) [`get_quadratic_quench`](@ref)

# Examples
```jldoctest
julia> thermostat = ProtoSyn.Drivers.get_constant_temperature(1.0)
(::ProtoSyn.Drivers.var"#constant_temperature#9"{Float64}) (generic function with 1 method)

julia> thermostat(100)
1.0
```
"""
function get_constant_temperature(init_temp::T) where {T <: AbstractFloat}
    return function constant_temperature(step::Int)
        return init_temp
    end
end

"""
    get_linear_quench(init_temp::T, n_steps::Int) where {T <: AbstractFloat}

Return a linear temperature quench thermostat function, beggining at temperature
`init_temp` and reducing linearly to 0.0 during `n_steps`.

# See also
[`get_constant_temperature`](@ref) [`get_quadratic_quench`](@ref)

# Examples
```jldoctest
julia> thermostat = ProtoSyn.Drivers.get_linear_quench(1.0, 100)
(::ProtoSyn.Drivers.var"#linear_quench#10"{Float64, Int64}) (generic function with 1 method)

julia> thermostat(0)
1.0

julia> thermostat(50)
0.5

julia> thermostat(100)
0.0
```
"""
function get_linear_quench(init_temp::T, n_steps::Int) where {T <: AbstractFloat}
    return function linear_quench(step::Int)
        return -((step * init_temp)/n_steps) + init_temp
    end
end

"""
    get_quadratic_quench(init_temp::T, n_steps::Int, min_temp::T) where {T <: AbstractFloat}

Return a quadratic temperature quench thermostat function, beginning at
temperature `init_temp` and reducing to `min_temp` during `n_steps` following a
quadratic curve.

# See also
[`get_constant_temperature`](@ref) [`get_linear_quench`](@ref)

# Examples
```jldoctest
julia> thermostat = ProtoSyn.Drivers.get_quadratic_quench(1.0, 100, 0.0)
(::ProtoSyn.Drivers.var"#quadratic_quench#15"{Float64, Int64, Float64}) (generic function with 1 method)

julia> thermostat(0)
1.0

julia> thermostat(50)
0.25

julia> thermostat(100)
0.0
```
"""
function get_quadratic_quench(init_temp::T, n_steps::Int, min_temp::T) where {T <: AbstractFloat}
    return function quadratic_quench(step::Int)
        a = (init_temp-min_temp)/(n_steps * n_steps)
        return a * (step - n_steps)^2 + min_temp
    end
end