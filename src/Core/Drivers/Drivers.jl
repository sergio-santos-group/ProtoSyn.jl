module Drivers

using ProtoSyn

abstract type DriverState end

abstract type Driver end

include("callback.jl")
include("steepest_descent.jl")
include("monte_carlo.jl")
include("ILS.jl")
include("compound.jl")
include("thermostats.jl")

# TODO - Molecular Dynamics
#include("md.jl")

end