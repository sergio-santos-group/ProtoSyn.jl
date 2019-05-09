module Forcefield

using ..Aux
using ..Common
using LinearAlgebra, JSON

include("Amber/Amber.jl")
include("Restraints/Restraints.jl")
include("CoarseGrain/CoarseGrain.jl")

evaluate!(s::Common.State, t::Amber.Topology, f::Bool = false) = Amber.evaluate!(s, t, f)
evaluate!(s::Common.State, t::Vector{Restraints.DistanceFBR}, f::Bool = false) = Restraints.evaluate!(s, t, f)

end