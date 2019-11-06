module ProtoSyn

using LinearAlgebra
#Statistics

const resource_dir = joinpath(dirname(@__DIR__), "resources")

#region CORE ------------------------------------
include("XMLRPC/XMLRPC.jl")

include("core/constants.jl")
include("core/macros.jl")
include("core/types.jl")
include("core/base.jl")
include("core/io.jl")
include("core/math.jl")

#endregion


#function eval! end
function run! end

#region SUBMODULES ------------------------------

include("Peptides/Peptides.jl")
# include("Forcefields/Forcefields.jl")
# include("Calculators/Calculators.jl")
# include("Drivers/Drivers.jl")

#endregion



include("nbdisplay.jl")
function __init__()
    if isdefined(Main, :IJulia) && Main.IJulia.inited
        println("Embedding javascript")
        embed_javascript()
    end
end


end # module

