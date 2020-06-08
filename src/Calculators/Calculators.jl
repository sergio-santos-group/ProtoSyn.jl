module Calculators

using ..ProtoSyn
using ..ProtoSyn: @dot, @cross
using ..ProtoSyn.Units: tonumber

# resource directory for this module
const resource_dir = let
    modname = string(nameof(@__MODULE__))
    joinpath(ProtoSyn.resource_dir, modname)
end


# hook point for function overloading
export eval!
function eval! end

include("state.jl")
include("types.jl")
include("forcefield.jl")

include("potentials.jl")
include("methods.jl")
include("bonded.jl")
include("nonbonded.jl")

#include("Forcefield/Forcefield.jl")

end
