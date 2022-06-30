@doc """
    Peptides

The Peptides modules introduces `Calculators`, `Mutators` and `Drivers` (among
other methods) specific for proteins and peptides.
"""
module Peptides

using ..ProtoSyn

# Resource directory for this module
const resource_dir = let
    modname = string(nameof(@__MODULE__))
    joinpath(ProtoSyn.resource_dir, modname)
end

include("Types/types.jl")

# Selections
include("Submodules/Selections/polar.jl")
include("Submodules/Selections/sidechains.jl")
include("Submodules/Selections/chi.jl")
include("Submodules/Selections/phi-psi-omega.jl")
include("Submodules/Selections/protein.jl")
include("Submodules/Selections/secondary-structure.jl") # Requires types.jl

include("constants.jl")
include("Submodules/Rotamers/Rotamers.jl")
include("Methods/io.jl")
include("Methods/graph.jl")
include("Methods/travel-graph.jl")
include("Methods/state.jl")
include("Methods/pose.jl")
include("Methods/other.jl") # Requires constants.jl
include("Submodules/Builder/Builder.jl")
include("Submodules/Builder/grammar.jl")
include("Calculators/Calculators.jl") # Requires grammar.jl
include("Mutators/Mutators.jl") # Requires Rotamers.jl
include("Drivers/drivers.jl")

# External packages
include("Submodules/ExternalPackages/GROMACS/gromacs.jl")

end