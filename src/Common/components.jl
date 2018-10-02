abstract type AbstractEnergy end

struct NullEnergy <: AbstractEnergy end


mutable struct State
    
    size::Int64
    energy::AbstractEnergy
    xyz::Array{Float64, 2}
    forces::Array{Float64, 2} # kJ mol⁻¹ nm⁻¹
    atnames::Array{String, 1}

end
State(n::Int64) = State(n, NullEnergy(), zeros(n, 3), zeros(n, 3), Array{String, 1}())


mutable struct Residue
    atoms::Array{Int64, 1}
    next::Union{Residue, Int64, Nothing}
    name::String
end


