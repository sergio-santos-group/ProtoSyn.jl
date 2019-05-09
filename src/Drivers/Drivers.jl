module Drivers

using ..Aux
using ..Common
using ..Abstract
using Printf

macro callback(freq::Int, ex::Expr)

    new_fcn_decl = :($(gensym())(
        a::Common.State,
        b::Abstract.DriverState,
        c::Abstract.DriverConfig)
    )

    fcn_decl = ex.args[1]
    if fcn_decl.head == :tuple
        # this is a nameless function:
        #  retain generated name
        start = 1
    else
        # this is a named function:
        # copy function name to new declaration
        new_fcn_decl.args[1] = fcn_decl.args[1]
        start = 2
    end
    
    varnames = map(s->s isa Expr ? s.args[1] : s, fcn_decl.args[start:end])
    varnames = varnames[1:min(3,length(varnames))]

    for (index,varname) in enumerate(varnames)
        new_fcn_decl.args[index+1].args[1] = varname
    end
    ex.args[1] = new_fcn_decl

    # ex
    :(Common.@callback $(freq) $(ex))
end

include("SteepestDescent/SteepestDescent.jl")
include("MonteCarlo/MonteCarlo.jl")
include("ILSRR/ILSRR.jl")
include("MD/MD.jl")

function Base.show(io::IO, b::Union{Abstract.DriverConfig, Abstract.DriverState})
    print(io, string(typeof(b)))
    for p in fieldnames(typeof(b))
        print(io, "\n   $(String(p)) = $(getproperty(b,p))")
    end
end

end