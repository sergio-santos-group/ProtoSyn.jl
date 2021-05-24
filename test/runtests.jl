push!(LOAD_PATH, "../src")

using Test, ProtoSyn, Printf

# include("graph.jl")
include("Core/units.jl")
include("Core/builder-submodule.jl")

# * The following code was previously tested on builder-submodule.jl and will be
# * used in the remaining test sets
res_lib = ProtoSyn.Peptides.grammar(Float64)
pose    = ProtoSyn.build(res_lib, seq"GME"); sync!(pose)
backup  = copy(pose)

include("Core/pose-methods.jl")
include("Core/state-methods.jl")
include("Core/graph-methods.jl")
include("Core/selections-submodule.jl")
include("Core/calculators.jl")
include("Core/mutators.jl")
include("Core/drivers.jl")

