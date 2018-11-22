module Common

using ..Aux

include("energy.jl")
include("residue.jl")
include("dihedral.jl")
include("metadata.jl")
include("state.jl")
include("callback.jl")
include("loaders.jl")
include("conf_gen.jl")
include("macros.jl")

const ss2bbd = Dict(
        Common.SS.SHEET => Dict(Common.DIHEDRAL.phi => deg2rad(-139.0), Common.DIHEDRAL.psi => deg2rad(135.0)),
        Common.SS.HELIX => Dict(Common.DIHEDRAL.phi => deg2rad(-57.0),  Common.DIHEDRAL.psi => deg2rad(-47.0)))
end