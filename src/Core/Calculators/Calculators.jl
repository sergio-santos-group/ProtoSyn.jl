module Calculators

using ProtoSyn
using Base.Cartesian

@time begin
    include("verlet_list.jl")
    include("distance_matrix.jl")

    # Load energy function components
    mutable struct EnergyFunctionComponent

        name::String
        calc::Function
        cached_n_calls::Int16
    end

    EnergyFunctionComponent(name::String, calc::Function) = begin
        return EnergyFunctionComponent(name, calc, Int16(0))
    end

    @info " | Loading TorchANI"
    include("torchani.jl")

    @info " | Loading Restraint Models"
    include("restraints.jl")

    @info " | Loading Energy Function"
    include("energy_function.jl")
end

end