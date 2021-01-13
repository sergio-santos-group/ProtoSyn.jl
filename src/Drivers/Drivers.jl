module Drivers

using ProtoSyn
using ProtoSyn: DriverState, Driver

include("callback.jl")
include("steepest_descent.jl")

# molecular dynamics
include("thermostats.jl")
#include("md.jl")

include("monte_carlo.jl")
include("compound.jl")
#(driver::IDriver)(s::State)  = driver(nothing, s)


end
