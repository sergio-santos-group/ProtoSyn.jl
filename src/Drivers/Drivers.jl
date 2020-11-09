module Drivers

using ProtoSyn

abstract type DriverState end
abstract type Driver end

#abstract type AbstractSampler end

# steepest descent
include("sd.jl")

# molecular dynamics
#include("thermostats.jl")
#include("md.jl")

include("mc2.jl")

#(driver::IDriver)(s::State)  = driver(nothing, s)


end
