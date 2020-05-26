module ProtoSyn

const resource_dir = joinpath(dirname(@__DIR__), "resources")

#region CORE ------------------------------------
include("XMLRPC/XMLRPC.jl")


# #endregion
include("Units.jl")
using .Units: tonumber

include("core/graph.jl")
include("core/macros.jl")
include("core/types.jl")
include("core/state.jl")

export Pose
struct Pose{T<:AbstractContainer}
    graph::T
    state::State
    Pose(c::T, s::State) where {T<:AbstractContainer}= begin
        c.id != s.id && error("unpairable container and state")
        new{T}(c,s)
    end
end
Base.copy(p::Pose) = Pose(copy(p.graph),copy(p.state))

include("core/base.jl")
include("core/selection.jl")


export Fragment
const Fragment = Pose{Segment}

export ResidueDB
const ResidueDB = Dict{String, Fragment}

include("core/io.jl")           # <-- ATTENTION
include("core/iterators.jl")    # <-- ATTENTION
include("core/methods.jl")
# include("core/loaders.jl")



#region SUBMODULES ------------------------------
include("Core/Builder.jl")

include("Peptides/Peptides.jl")
include("Sugars/Sugars.jl")
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

