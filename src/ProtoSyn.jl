module ProtoSyn

__precompile__(true)

const resource_dir = joinpath(dirname(@__DIR__), "resources")

@info "Loading Core"
@time begin
    include("Core/XMLRPC/XMLRPC.jl")
    include("Core/Units/Units.jl")
    include("Core/Methods/macros.jl")
    include("Core/Types/common.jl")
    include("Core/Types/graph.jl")
    include("Core/Types/state.jl")
    include("Core/Types/pose.jl")
    include("Core/Methods/graph.jl")
    include("Core/Methods/state.jl")
    include("Core/Methods/base.jl")
    include("Core/Methods/io.jl")
    include("Core/Methods/iterators.jl")
    include("Core/Selections/selections.jl")
    include("Core/Builder/grammar.jl")
    include("Core/Builder/Builder.jl")
end
Base.copy(p::Pose) = Pose(copy(p.graph),copy(p.state))

include("Core/base.jl")
include("Core/selection2.jl")


export Fragment
const Fragment = Pose{Segment}

export ResidueDB
const ResidueDB = Dict{String, Fragment}

include("Core/io.jl")           # <-- ATTENTION
include("Core/iterators.jl")    # <-- ATTENTION
include("Core/methods.jl")
# include("Core/loaders.jl")



#region SUBMODULES ------------------------------
include("Core/Builder/Builder.jl")

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

@info "Loading Peptides"
@time include("Peptides/Peptides.jl")

# include("Sugars/Sugars.jl")

@info "Loading Calculators"
include("Calculators/Calculators.jl")

@info "ProtoSyn loaded successfully!"
end # module

