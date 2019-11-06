# ----------------------------------------------------------------------------------------------------------
#                                                 ENERGY
@doc raw"""
    Energy(eTotal::Float64)

Simple energy container.

# Examples
```julia-repl
julia> Common.Energy(1.15)
Energy(eTotal=1.15)

julia> Common.Energy()
Energy(eTotal=0.0)
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




function set_energy_component(container::Energy, comp::Symbol, value::Float64)
    container.components[comp] = value
    container.total = value
end
