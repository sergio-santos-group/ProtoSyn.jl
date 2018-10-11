module ProtoSyn

using LinearAlgebra, JSON, Printf

export Common, Aux, Forcefield, Print, Mutators, Drivers

include("Common/Common.jl")
include("Aux/Aux.jl")
include("Forcefield/Forcefield.jl")
include("Print/Print.jl")
include("Drivers/Drivers.jl")   # Order ?
include("Mutators/Mutators.jl") # Order ?

end