module Forcefield

using ...ProtoSyn
using ..Calculators
using ...ProtoSyn.Units: tonumber
using ...ProtoSyn: @dot, @cross

# resource directory for this module
const resource_dir = let
    modname = string(nameof(@__MODULE__))
    joinpath(Calculators.resource_dir, modname)
end



include("state.jl")
include("types.jl")
include("potentials.jl")
include("methods.jl")
include("bonded.jl")
include("nonbonded.jl")


# include("macros.jl")
# include("nonbonded.jl")



end