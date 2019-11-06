module Calculators

using ..ProtoSyn

# resource directory for this module
const resource_dir = let
    modname = string(nameof(@__MODULE__))
    joinpath(ProtoSyn.resource_dir, modname)
end


# hook point for function overloading
function eval! end

include("Forcefield/Forcefield.jl")

end
