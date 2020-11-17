module Drivers

using ProtoSyn

abstract type DriverState end
abstract type Driver end

#abstract type AbstractSampler end

include("callback.jl")
include("blitz.jl")

# steepest descent
include("sd.jl")

# molecular dynamics
#include("thermostats.jl")
#include("md.jl")

include("monte_carlo.jl")
include("compound.jl")
#(driver::IDriver)(s::State)  = driver(nothing, s)


end
