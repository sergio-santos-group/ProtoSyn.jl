module Calculators

using ProtoSyn
using Base.Cartesian

include("verlet_list.jl")
include("distance_matrix.jl")

# Load energy function components
mutable struct EnergyFunctionComponent

    name::String
    calc::Function
end

@info " | Loading TorchANI"
include("torchani.jl")

@info " | Loading Restraint Models"
include("restraints.jl")

@info " | Loading Energy Function"
include("energy_function.jl")

end