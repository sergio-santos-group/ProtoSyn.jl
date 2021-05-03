push!(LOAD_PATH, "../src")

using Base.Test, ProtoSyn

# include("graph.jl")
include("Core/units.jl")
include("Core/builder-submodule.jl")
include("Core/pose-methods.jl")
include("Core/state-methods.jl")
include("Core/graph-methods.jl")
include("Core/selections-submodule.jl")
include("Core/calculators.jl")

