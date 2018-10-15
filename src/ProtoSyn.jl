module ProtoSyn

using LinearAlgebra, JSON, Printf

export Common, Aux, Forcefield, Print, Mutators, Drivers

include("Aux/Aux.jl")
include("Common/Common.jl")
include("Forcefield/Forcefield.jl")
include("Print/Print.jl")
include("Mutators/Mutators.jl") # Order ?
include("Drivers/Drivers.jl")   # Order ?

end