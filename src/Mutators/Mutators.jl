module Mutators

using ..Forcefield
using ..Drivers
using ..Common
using ..Aux
using ..Print
using ..Abstract

# abstract type AbstractMutator end

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

apply!(st::Common.State, mut::Dihedral.MutatorConfig) = Dihedral.apply!(st, mut)
# apply!(st::Common.State, mut::Crankshaft.CrankshaftMutator) = Crankshaft.apply!(st, mut)

mutable struct Sampler{F <: Function, G <: Function, T <: Abstract.MutatorConfig, N}

    # Parameters:            Signatures:
    apply!::F                # sampler.apply!(state::Common.State, sampler.mutators::Vector{Abstract.MutatorConfig})
    tune!::Union{G, Nothing} # sampler.tune!(sampler.mutators::Vector{Abstract.MutatorConfig}, driver_state::AbstractDriverState)
    mutators::Tuple          # Note: Can't add type, but type specificiy is assured by the constructor
end # mutable struct

# -----------------------------------
# Option A:
# function aggregate!(state::Common.State, mutators::Tuple{T}) where {T <: Abstract.MutatorConfig}
#     for mutator in mutators
#         Mutators.apply!(state, mutator)
#     end # for
# end # function

# function Sampler(mutators::Vararg{T}) where {T <: Abstract.MutatorConfig}
#     Sampler{Function, Function}(aggregate!, nothing, mutators)
# end

# function Sampler(apply!::F, mutators::Vararg{T}) where {F <: Function, T <: Abstract.MutatorConfig}
#     Sampler{Function, Function}(apply!, nothing, mutators)
# end

# function Sampler(apply!::F, tune!::G, mutators::Vararg{T}) where {F <: Function, G <: Function, T <: Abstract.MutatorConfig}
#     Sampler{Function, Function}(apply!, tune!, mutators)
# end

# # Exameple calls:
# my_sampler! = Mutators.Sampler(rand, randn, dihedral_mutator, dihedral_mutator)
# my_sampler! = Mutators.Sampler(rand, dihedral_mutator)
# my_sampler! = Mutators.Sampler(dihedral_mutator); my_sampler!.tune! = rand

# -----------------------------------
# Option B:
function Sampler(mutators::T...; apply!::Union{F, Nothing} = nothing, tune!::Union{G, Nothing} = nothing) where { T <: Abstract.MutatorConfig, F <: Function, G <: Function}
    if apply! == nothing
        apply! = function aggregate!(state::Common.State, mutators::T...) where {T <: Abstract.MutatorConfig}
            for mutator in mutators
                Mutators.apply!(state, mutator)
            end
        end
    end
    Sampler{Function, Function, Abstract.MutatorConfig, Int64}(apply!, tune!, mutators)
end

# Example calls:
# my_sampler! = Mutators.Sampler(dihedral_mutator, dihedral_mutator, tune! = rand, apply! = rand)
# my_sampler! = Mutators.Sampler(tune! = randn, apply! = rand)
# my_sampler! = Mutators.Sampler(dihedral_mutator, tune! = rand)


function Base.show(io::IO, b::Sampler)
    tuner = b.tune! == nothing ? "Nothing" : string(b.tune!)
    print(io, "Sampler(\n mutators=$(length(b.mutators)),\n tune!=$(tuner),\n apply!=$(string(b.apply!))\n)")
end

end # module