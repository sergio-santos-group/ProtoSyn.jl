module ProtoSyn

const resource_dir = joinpath(dirname(@__DIR__), "resources")

#region CORE ------------------------------------
include("Core/XMLRPC/XMLRPC.jl")
include("Core/Units/Units.jl")
using .Units: tonumber


# #endregion

include("Core/graph.jl")
include("Core/macros.jl")
include("Core/types.jl")
include("Core/state.jl")

export Pose
struct Pose{T<:AbstractContainer}
    graph::T
    state::State
    Pose(c::T, s::State) where {T<:AbstractContainer}= begin
        c.id != s.id && error("unpairable container and state")
        new{T}(c, s)
    end
end
Base.copy(p::Pose) = Pose(copy(p.graph),copy(p.state))

include("Core/base.jl")


export Fragment
const Fragment = Pose{Segment}

export ResidueDB
const ResidueDB = Dict{String, Fragment}

include("Core/io.jl")           # <-- ATTENTION
include("Core/iterators.jl")    # <-- ATTENTION
include("Core/methods.jl")

include("Core/Selections/selections.jl") # Makes use of iterators, must come after

# include("Core/loaders.jl")



#region SUBMODULES ------------------------------
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

