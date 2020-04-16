module Drivers

using ..ProtoSyn

# IDriverState
#   defines:
#   - step::Int
#   - stalled::Bool
#   - converged::Bool
#   - completed::Bool
abstract type IDriverState end


# IDriver
#   defines:
#   - eval!::Function
#   - max_steps::Int
#   - pairlist_freq::Int
abstract type IDriver end

#abstract type AbstractSampler end

# steepest descent
#include("sd.jl")

# molecular dynamics
#include("thermostats.jl")
#include("md.jl")

include("mc2.jl")

#(driver::IDriver)(s::State)  = driver(nothing, s)


end
