module Mutators

using ..Forcefield
using ..Drivers
using ..Common
using ..Aux
using ..Print

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

end