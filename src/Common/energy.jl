# ----------------------------------------------------------------------------------------------------------
#                                                 ENERGY
@doc raw"""
    Energy(; components::Dict{Symbol, Float64} = Dict{Symbol, Float64}(), total::Float64 = 0.0)

Simple energy container.

# Examples
```julia-repl
julia> Common.Energy(total = 1.15)
Energy(eTotal=1.15)

julia> Common.Energy()
Energy(total=0.0, components=Dict{Symbol,Float64}())
```
"""
Base.@kwdef mutable struct Energy
    components::Dict{Symbol, Float64} = Dict{Symbol, Float64}()
    total::Float64             = 0.0
end
Base.show(io::IO, b::Energy) = print(io, "Energy(total=$(b.total), components=$(b.components))")

function Base.copy!(dst::Energy, src::Energy)::Energy
    copy!(dst.components, src.components)
    dst.total = src.total
    return dst
end


@doc raw"""
    set_energy_component!(container::Energy, comp::Symbol, value::Float64)

Updates the component in Energy `container`; updates Energy.total to reflect the given `value`.

# Examples
```julia-repl
julia> Common.set_energy_component!(state.energy, :bonds, 10.5)
```
"""
function set_energy_component!(container::Energy, comp::Symbol, value::Float64)
    container.components[comp] = value
    container.total = value
end