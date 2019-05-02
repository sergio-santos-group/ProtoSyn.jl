module Drivers

# using ..Aux
# using ..Common
# using Printf


abstract type AbstractDriverConfig end
abstract type AbstractDriverState  end

include("SteepestDescent/SteepestDescent.jl")
include("MonteCarlo/MonteCarlo.jl")
include("ILSRR/ILSRR.jl")
include("MD/MD.jl")

function Base.show(io::IO, b::Union{AbstractDriverConfig, AbstractDriverState})
    print(io, string(typeof(b)))
    for p in fieldnames(typeof(b))
        print(io, "\n   $(String(p)) = $(getproperty(b,p))")
    end
end

end