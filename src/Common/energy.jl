# ----------------------------------------------------------------------------------------------------------
#                                                 ENERGY

abstract type AbstractEnergy end

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
Base.@kwdef mutable struct Energy <: AbstractEnergy
    comp::Dict{String, Float64} = Dict{String, Float64}()
    eTotal::Float64             = 0.0
end
Base.show(io::IO, b::Energy) = print(io, "Energy(eTotal=$(b.eTotal), components=$(b.comp))")

function Base.copy!(dst::Energy, src::Energy)::Energy
    copy!(dst.comp, src.comp)
    dst.eTotal = src.eTotal
    return dst
end