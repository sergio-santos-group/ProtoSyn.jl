module Calculators

using ProtoSyn
using Base.Cartesian
using Printf

include("verlet_list.jl")
include("distance_matrix.jl")

# Load energy function components
mutable struct EnergyFunctionComponent{T <: AbstractFloat}

    name::String
    calc::Function
    settings::Dict{Symbol, Any}
    α::T
    update_forces::Bool
end

function Base.show(io::IO, efc::EnergyFunctionComponent{T}) where {T <: AbstractFloat}
    @printf(io, "%14s : %-s\n", "Name", efc.name)
    @printf(io, "%14s : %-s\n", "Weight (α)", efc.α)
    @printf(io, "%14s : %-s\n", "Update forces", string(efc.update_forces))
    if length(efc.settings) == 0
        @printf(io, "%14s : -", "Setings")
    else
        @printf(io, "%14s :\n", "Setings")
        for (key, value) in efc.settings
            @printf(io, "%15s => ", ":"*string(key))
            if typeof(value) == Matrix{T}
                print(io, "Matrix{$T}($(size(value))\n")
            else
                print(io, "$value\n")
            end
        end
    end
end

@info " | Loading TorchANI"
include("torchani.jl")

@info " | Loading Restraint Models"
include("potential.jl")
include("restraints.jl")

@info " | Loading Energy Function"
include("energy_function.jl")

end