module Mutators

using ..Forcefield
using ..Drivers
using ..Common
using ..Aux
using ..Print
using ..Abstract

include("Dihedral/Dihedral.jl")
include("Crankshaft/Crankshaft.jl")
include("Blockrot/Blockrot.jl")
include("Sidechain/Sidechain.jl")

apply!(st::Common.State, mut::Dihedral.MutatorConfig)   = Dihedral.apply!(st, mut)
apply!(st::Common.State, mut::Crankshaft.MutatorConfig) = Crankshaft.apply!(st, mut)
apply!(st::Common.State, mut::Blockrot.MutatorConfig)   = Blockrot.apply!(st, mut)
apply!(st::Common.State, mut::Sidechain.MutatorConfig)  = Sidechain.apply!(st, mut)

apply!(st::Common.State, mut::Drivers.SteepestDescent.DriverConfig) = Drivers.SteepestDescent.run!(st, mut)

@doc raw"""
    Sampler(; mutators::Vector{Any} = [], apply! = Union{Function, Nothing} = nothing, tune! = Union{Function, Nothing} = nothing) <: Abstract.Sampler

A Sampler is an aggregator of Abstract.MutatorConfigs and/or Abstract.DriverConfigs (`mutators`), who are applied according to an `apply!` function.
If no function is passed to the constructor, the default function is simply a for loop running over all mutators, in order.
A `tune!` may be passed, which is called in certain drivers and operates over the parameters of the `mutators` to tune them in
response to changes during the runtime of the Driver. 

# Examples
```julia-repl
julia> Forcefield.Sampler(mutators = [dihedral, crankshaft])
```
"""
mutable struct Sampler{F <: Function, G <: Function, T <: Any} <: Abstract.Sampler

    # Parameters:            Signatures:
    apply!::F                # sampler.apply!(state::Common.State, sampler.mutators::Vector{Abstract.MutatorConfig})
    tune!::Union{G, Nothing} # sampler.tune!(sampler.mutators::Vector{Abstract.MutatorConfig}, driver_state::Abstract.DriverState)
    mutators::Vector{T}
end # end struct

function Sampler(; mutators::Vector{T} = Vector{T}(), apply!::Union{F, Nothing} = nothing, tune!::Union{G, Nothing} = nothing) where {F <: Function, G <: Function, T <: Any}
    if apply! == nothing
        apply! = function default_aggregate!(state::Common.State, mutators::Vector{T}) where {T <: Any}
            for mutator in mutators
                Mutators.apply!(state, mutator)
            end # end for
        end # end function
    end # end if
    Sampler{Function, Function, Any}(apply!, tune!, mutators)
end # end function

function Base.show(io::IO, b::Sampler)
    tuner = b.tune! == nothing ? "Nothing" : string(b.tune!)
    print(io, "Sampler(\n mutators=$(length(b.mutators)),\n tune!=$(tuner),\n apply!=$(string(b.apply!))\n)")
end # end function

end # module