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
mutable struct Energy <: AbstractEnergy
    eTotal::Float64
end
Energy() = Energy(0.0)
Base.show(io::IO, b::Energy) = print(io, "Energy(eTotal=$(b.eTotal))")