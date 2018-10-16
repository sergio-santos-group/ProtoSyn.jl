module ProtoSyn

using LinearAlgebra, JSON, Printf

export Common, Aux, Forcefield, Print, Mutators, Drivers

include("Common/Common.jl")
include("Aux/Aux.jl")
include("Forcefield/Forcefield.jl")
include("Print/Print.jl")
include("Mutators/Mutators.jl")
include("Drivers/Drivers.jl")

include("nbdisplay.jl")

function __init__()
    if isdefined(Main, :IJulia) && Main.IJulia.inited
        embed_javascript()
    end
end

end