module Forcefield

using ..Aux
using ..Common
using LinearAlgebra, JSON

include("Amber/Amber.jl")
include("Restraints/Restraints.jl")
include("CoarseGrain/CoarseGrain.jl")

evaluate!(t::Amber.Topology, s::Common.State; cut::Float64 = 2.0, do_forces = false) = Amber.evaluate!(t,s,cut,do_forces)
evaluate!(t::Vector{Restraints.DistanceFBR}, s::Common.State; do_forces::Bool = false) = Restraints.evaluate!(t,s,do_forces)

end