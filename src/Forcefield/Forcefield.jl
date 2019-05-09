module Forcefield

using ..Aux
using ..Common
using ..Abstract
using LinearAlgebra, JSON

include("Amber/Amber.jl")
include("Restraints/Restraints.jl")
include("CoarseGrain/CoarseGrain.jl")

evaluate!(s::Common.State, t::Vector{Amber.HarmonicBond},     f::Bool = false) = Amber.evaluate!(s, t, f)
evaluate!(s::Common.State, t::Vector{Amber.HarmonicAngle},    f::Bool = false) = Amber.evaluate!(s, t, f)
evaluate!(s::Common.State, t::Vector{Amber.DihedralCos},      f::Bool = false) = Amber.evaluate!(s, t, f)
evaluate!(s::Common.State, t::Vector{Amber.Atom},             f::Bool = false) = Amber.evaluate!(s, t, f)
evaluate!(s::Common.State, t::Amber.Topology,                 f::Bool = false) = Amber.evaluate!(s, t, f)

evaluate!(s::Common.State, t::Vector{Restraints.DistanceFBR}, f::Bool = false) = Restraints.evaluate!(s, t, f)
evaluate!(s::Common.State, t::Vector{Restraints.DihedralFBR}, f::Bool = false) = Restraints.evaluate!(s, t, f)

evaluate!(s::Common.State, t::Vector{CoarseGrain.SolvPair},   f::Bool = false) = CoarseGrain.evaluate!(s, t, f)
evaluate!(s::Common.State, t::Vector{CoarseGrain.HbNetwork},  f::Bool = false) = CoarseGrain.evaluate!(s, t, f)


# TODO: Documentation
mutable struct Evaluator{F <: Function, T <: Abstract.ForcefieldComponent} <: Abstract.Evaluator

    # Parameters:            Signatures:
    evaluate!::F             # evaluator.evaluate!(state::Common.State, evaluator.components::Vector{Abstract.ForcefieldComponent}, do_forces::Bool = false)
    components::Vector{T}
end # mutable struct


# TODO: Documentation
function Evaluator(; components::Vector{T} = Vector{T}(), evaluate!::Union{F, Nothing} = nothing) where {F <: Function, T <: Abstract.ForcefieldComponent}
    if evaluate! == nothing
        evaluate! = function default_sum!(state::Common.State, components::Vector{T}, do_forces::Bool = false)::Float64 where {T <: Abstract.ForcefieldComponent}
            total = 0.0
            for component in components
                total += Forcefield.evaluate!(state, component, do_forces)
            end
            state.energy.total = total
            return total
        end
    end
    Evaluator{Function, Abstract.ForcefieldComponent}(evaluate!, components)
end

end