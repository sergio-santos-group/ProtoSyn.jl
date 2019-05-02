# ----------------------------------------------------------------------------------------------------------
#                                                 ENERGY

# abstract type AbstractEnergy end
abstract type AbstractEnergy{T<:AbstractFloat} end

@doc raw"""
    NullEnergy()

Empty placeholder energy container.

# Examples
```julia-repl
julia> Common.NullEnergy()
Null
```
"""
struct NullEnergy <: AbstractEnergy end
Base.show(io::IO, b::NullEnergy) = print(io, "Null")


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
mutable struct Energy <: AbstractEnergy
    comp::Dict{String, Float64}
    eTotal::Float64
end
Energy() = Energy(Dict{String, Float64}("other" => 0.0, "amber" => 0.0), 0.0)
Base.show(io::IO, b::Energy) = print(io, "Energy(eTotal=$(b.eTotal), components=$(b.comp))")

function Base.copy!(dst::Energy, src::Energy)::Energy
    copy!(dst.comp, src.comp)
    dst.eTotal = src.eTotal
    return dst
end



function set_energy_component(container::AbstractEnergy, comp::Symbol, value::Float64)
    setproperty!(container, comp, value)
    container.total = value
end

