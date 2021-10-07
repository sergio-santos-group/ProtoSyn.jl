module Calculators

using ProtoSyn
using Base.Cartesian
using Printf

include("verlet_list.jl")
include("distance_matrix.jl")

# Load energy function components
include("energy_function_component.jl")

if "NO_TORCHANI" in keys(ENV) & ENV["NO_TORCHANI"] === "true"
    @warn "Environment variable NO_TORCHANI set to `true`. Not loading torchani."
else
    @info " | Loading TorchANI"
    include("torchani.jl")
end

@info " | Loading Restraint Models"
include("Potentials/potentials.jl")
include("restraints.jl")

@info " | Loading Energy Function"
include("energy_function.jl")

end