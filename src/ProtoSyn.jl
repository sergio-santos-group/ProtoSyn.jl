module ProtoSyn

const resource_dir = joinpath(dirname(@__DIR__), "resources")

#region CORE ------------------------------------
include("Core/XMLRPC/XMLRPC.jl")
include("Core/Units/Units.jl")
using .Units: tonumber


# #endregion CORE

include("Core/Methods/macros.jl")
include("Core/Types/graph.jl")
include("Core/Types/state.jl")
include("Core/Types/pose.jl")
include("Core/Methods/graph.jl")
include("Core/Methods/state.jl")

include("Core/Methods/base.jl")



# Fragment(pose::Pose) = begin
#     Fragment()
# end

export ResidueDB
const ResidueDB = Dict{String, Fragment}

include("Core/Methods/io.jl")           # <-- ATTENTION
include("Core/Methods/iterators.jl")    # <-- ATTENTION

include("Core/Selections/selections.jl") # Makes use of iterators, must come after

# include("Core/loaders.jl")



#region SUBMODULES ------------------------------
include("Core/Builder/grammar.jl")
include("Core/Builder/Builder.jl")
include("Peptides/Peptides.jl")

include("Sugars/Sugars.jl")
# include("Forcefields/Forcefields.jl")
include("Calculators/Calculators.jl")
include("Drivers/Drivers.jl")

#endregion



# include("nbdisplay.jl")
# function __init__()
#     if isdefined(Main, :IJulia) && Main.IJulia.inited
#         println("Embedding javascript")
#         embed_javascript()
#     end
# end




end # module

