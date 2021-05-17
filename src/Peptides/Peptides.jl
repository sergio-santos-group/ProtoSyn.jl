@doc """
    Peptides

The Peptides modules introduces `Calculators`, `Mutators` and `Drivers` (among
other methods) specific for proteins and peptides.
"""
module Peptides

using ..ProtoSyn

# resource directory for this module
const resource_dir = let
    modname = string(nameof(@__MODULE__))
    joinpath(ProtoSyn.resource_dir, modname)
end

# include("constants.jl")
# include("Calculators/Calculators.jl")
# include("Rotamers/Rotamers.jl")
# include("io.jl")
include("methods.jl")
# include("Mutators/Mutators.jl")
# include("Drivers/drivers.jl")

end