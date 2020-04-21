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


include("core/graph.jl")
include("core/macros.jl")
include("core/types.jl")
include("core/base.jl")
include("core/state.jl")

include("core/selection2.jl")

# export CSPair
# const CSPair = Tuple{AbstractContainer, State}
# pair(c::AbstractContainer, s::State) = begin
#     p = c
#     while hascontainer(p)
#         p = p.container
#     end
#     p.id != s.id && error("unable to pair container with state having different IDs")
#     Pair(c,s)
# end
# 
# @inline container(csp::CSPair) = csp[1]
# @inline state(csp::CSPair) = csp[2]
# 
# container(t::Topology; topmost::Bool) = t
# container(c::AbstractContainer; topmost=false) = begin
#     topmost && hascontainer(c) ? container(c.container; topmost=true) : c.container
# end


export ResidueDB
const ResidueDB = Dict{String, Tuple{Residue, State}}

include("core/io.jl")           # <-- ATTENTION
include("core/iterators.jl")    # <-- ATTENTION
include("core/methods.jl")


export from
from(db::ResidueDB, key::String) = begin
    r,s=db[key]
    copy(r), copy(s)
end

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

