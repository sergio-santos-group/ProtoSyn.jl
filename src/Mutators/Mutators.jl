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

# function Base.show(io::IO, b::Union{AbstractDriverConfig, AbstractDriverState})
#     print(io, string(typeof(b)))
#     for p in fieldnames(typeof(b))
#         print(io, "\n   $(String(p)) = $(getproperty(b,p))")
#     end
# end

apply!(st::Common.State, mut::Dihedral.MutatorConfig)   = Dihedral.apply!(st, mut)
apply!(st::Common.State, mut::Crankshaft.MutatorConfig) = Crankshaft.apply!(st, mut)
apply!(st::Common.State, mut::Blockrot.MutatorConfig)   = Blockrot.apply!(st, mut)
apply!(st::Common.State, mut::Sidechain.MutatorConfig)  = Sidechain.apply!(st, mut)

apply!(st::Common.State, mut::Drivers.SteepestDescent.DriverConfig) = Drivers.SteepestDescent.run!(st, mut)

# TODO: Documentation
mutable struct Sampler{F <: Function, G <: Function} <: Abstract.Sampler

    # Parameters:            Signatures:
    apply!::F                # sampler.apply!(state::Common.State, sampler.mutators::Vector{Abstract.MutatorConfig})
    tune!::Union{G, Nothing} # sampler.tune!(sampler.mutators::Vector{Abstract.MutatorConfig}, driver_state::Abstract.DriverState)
    mutators::Vector{Any}
end # mutable struct


# TODO Documentation
function Sampler(; mutators::Vector{Any} = Vector{T}(), apply!::Union{F, Nothing} = nothing, tune!::Union{G, Nothing} = nothing) where {F <: Function, G <: Function}
    if apply! == nothing
        apply! = function default_aggregate!(state::Common.State, mutators::Vector{Any})
            for mutator in mutators
                Mutators.apply!(state, mutator)
            end
        end
    end
    Sampler{Function, Function}(apply!, tune!, mutators)
end


function Base.show(io::IO, b::Sampler)
    tuner = b.tune! == nothing ? "Nothing" : string(b.tune!)
    print(io, "Sampler(\n mutators=$(length(b.mutators)),\n tune!=$(tuner),\n apply!=$(string(b.apply!))\n)")
end

end # module