module ProtoSyn

using LinearAlgebra, JSON, Printf, StatsBase

export Common, Aux, Forcefield, Print, Mutators, Drivers

include("Aux/Aux.jl")
include("Common/Common.jl")
include("Forcefield/Forcefield.jl")
include("Print/Print.jl")
include("Drivers/Drivers.jl")
include("Mutators/Mutators.jl")

include("nbdisplay.jl")

function __init__()
    if isdefined(Main, :IJulia) && Main.IJulia.inited
        embed_javascript()
    end
end

end