module ProtoSyn

# using LinearAlgebra
# #Statistics

const resource_dir = joinpath(dirname(@__DIR__), "resources")

# #region CORE ------------------------------------
include("XMLRPC/XMLRPC.jl")

# include("core/constants.jl")
# include("core/macros.jl")
# include("core/types.jl")
# include("core/base.jl")
# include("core/io.jl")
# include("core/math.jl")

# #endregion


# #function eval! end
# function run! end


include("core/macros.jl")
include("core/types.jl")
include("core/graph.jl")
include("core/base.jl")
include("core/state.jl")
include("core/io.jl")
include("core/iterators.jl")
include("core/methods.jl")

export ResidueDB
const ResidueDB = Dict{String, Tuple{Residue, State}}


#include("core/state2.jl")
#include("core/base2.jl")
#include("core/io2.jl")
#include("core/math2.jl")

#region SUBMODULES ------------------------------

include("Peptides/Peptides.jl")
# include("Forcefields/Forcefields.jl")
# include("Calculators/Calculators.jl")
# include("Drivers/Drivers.jl")

#endregion



# include("nbdisplay.jl")
# function __init__()
#     if isdefined(Main, :IJulia) && Main.IJulia.inited
#         println("Embedding javascript")
#         embed_javascript()
#     end
# end





end # module

