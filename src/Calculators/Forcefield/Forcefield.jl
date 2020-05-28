module Forcefield

# using ...ProtoSyn
# using ..Calculators
# using YAML

# resource directory for this module
# const resource_dir = let
#     modname = string(nameof(@__MODULE__))
#     joinpath(Calculators.resource_dir, modname)
# end


# include("types.jl")
include("methods.jl")
# include("macros.jl")

# include("bonded.jl")
# include("nonbonded.jl")



end