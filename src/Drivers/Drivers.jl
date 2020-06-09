module Drivers

using ..ProtoSyn

# IDriverState
#   defines:
#   - step::Int
#   - stalled::Bool
#   - converged::Bool
#   - completed::Bool
abstract type IDriverState end

# @inline step(s::IDriverState) = (s.step += 1; s)

# @inline stall(s::IDriverState) = (s.stalled = true; s)
# @inline complete(s::IDriverState) = (s.completed = true; s)
# @inline converge(s::IDriverState) = (s.converged = true; s)

# @inline isstalled(s::IDriverState) = s.stalled
# @inline isconverged(s::IDriverState) = s.converged
# @inline iscompleted(s::IDriverState) = s.completed


# IDriver
#   defines:
#   - eval!::Function
#   - max_steps::Int
#   - pairlist_freq::Int
abstract type IDriver end

#abstract type AbstractSampler end

# steepest descent
include("sd.jl")

# molecular dynamics
#include("thermostats.jl")
#include("md.jl")

include("mc2.jl")

#(driver::IDriver)(s::State)  = driver(nothing, s)


end
