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
evaluate!(s::Common.State, t::CoarseGrain.HbNetwork,          f::Bool = false) = CoarseGrain.evaluate!(s, t, f)


@doc raw"""
    Evaluator(; components::Vector{Any} = [], evaluate! = Union{Function, Nothing} = nothing) <: Abstract.Evaluator

An Evaluator is an aggregator of Abstract.ForcefieldComponents (`components`), who are evaluated according to an `evaluate!` function.
If no function is passed to the constructor, the default function is simply a sum of all components, in order.

# Examples
```julia-repl
julia> Forcefield.Evaluator(components = [amber_top, restraints])
```
"""
mutable struct Evaluator{F <: Function, T <: Any} <: Abstract.Evaluator

    # Parameters:            Signatures:
    evaluate!::F             # evaluator.evaluate!(state::Common.State, evaluator.components::Vector{Abstract.ForcefieldComponent}, do_forces::Bool = false)
    components::Vector{T}
end # end struct

function Evaluator(; components::Vector{T} = Vector{T}(), evaluate!::Union{F, Nothing} = nothing) where {F <: Function, T <: Any}
    if evaluate! == nothing
        evaluate! = function default_sum!(state::Common.State, components::Vector{T}, do_forces::Bool = false)::Float64 where {T <: Any}
            total = 0.0
            for component in components
                total += Forcefield.evaluate!(state, component, do_forces)
            end # end for
            state.energy.total = total
            return total
        end # end function
    end # end if
    Evaluator{Function, Any}(evaluate!, components)
end # end function

end # end module