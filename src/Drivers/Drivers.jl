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
    end # end if
    
    varnames = map(s->s isa Expr ? s.args[1] : s, fcn_decl.args[start:end])
    varnames = varnames[1:min(3,length(varnames))]

    for (index,varname) in enumerate(varnames)
        new_fcn_decl.args[index+1].args[1] = varname
    end # end for
    ex.args[1] = new_fcn_decl

    # ex
    :(Common.@callback $(freq) $(ex))
end # end macro

include("SteepestDescent/SteepestDescent.jl")
include("MonteCarlo/MonteCarlo.jl")
include("ILSRR/ILSRR.jl")
include("MolecularDynamics/MolecularDynamics.jl")

run!(s::Common.State, d::MonteCarlo.DriverConfig)      = MonteCarlo.run!(s, d)
run!(s::Common.State, d::SteepestDescent.DriverConfig) = SteepestDescent.run!(s, d)
run!(s::Common.State, d::ILSRR.DriverConfig)           = ILSRR.run!(s, d)

function Base.show(io::IO, b::Union{Abstract.DriverConfig, Abstract.DriverState})
    print(io, string(typeof(b)))
    for p in fieldnames(typeof(b))
        _p = getproperty(b,p) == nothing ? "nothing" : getproperty(b,p)
        print(io, "\n   $(String(p)) = $_p")
    end # end for
end # end function

end # end module