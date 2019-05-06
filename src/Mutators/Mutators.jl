module Mutators

using ..Forcefield
using ..Drivers
using ..Common
using ..Aux
using ..Print

abstract type AbstractMutator end

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

apply!(st::Common.State, mut::Dihedral.DihedralMutator) = Dihedral.apply!(st, mut)
apply!(st::Common.State, mut::Crankshaft.CrankshaftMutator) = Crankshaft.apply!(st, mut)

mutable struct Sampler{F <: Function} <: Function

    mutators::Vector{AbstractMutator}
    apply!::F

    Sampler(mutators::Any) = begin
        if all(isa.(mutators, AbstractMutator))
            new{Function}(mutators, function _apply!(state)
                for mutator in mutators
                    apply!(state, mutator)
                end #for
            end) # new
        else
            error("All mutators in Sampler constructor must be of type AbstractMutator")
        end # if
    end # Sampler
end # mutable struct
(sampler!::Sampler)(state) = sampler!.apply!(state)
Base.show(io::IO, b::Sampler) = print(io, "Sampler(mutators=$(length(b.mutators)))")

end # module