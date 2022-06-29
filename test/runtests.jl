push!(LOAD_PATH, "../src")

using Test, ProtoSyn, Printf
ProtoSyn.set_logger_to_error()

@testset verbose = true "Tests" begin

    # ------- CORE -------------------------------------------------------------
    include("Core/units.jl")
    include("Core/builder-submodule.jl")

    # * The following code was previously tested on builder-submodule.jl and
    # * will be used in the remaining test sets
    global res_lib = ProtoSyn.Peptides.grammar
    global pose    = ProtoSyn.build(res_lib, seq"GME"); sync!(pose)
    global backup  = copy(pose)

    include("Core/pose-methods.jl")
    include("Core/state-methods.jl")
    include("Core/graph-methods.jl")
    include("Core/selections-submodule.jl")
    include("Core/calculators.jl")
    include("Core/mutators.jl")
    include("Core/drivers.jl")

    # ------- PEPTIDES ---------------------------------------------------------
    include("Peptides/pose-methods.jl")
    include("Peptides/selections.jl")
    include("Peptides/rotamers.jl")

    # * The following code was previously tested on Peptides/rotamers.jl and
    # * will be used in the remaining test sets
    global rot_lib = ProtoSyn.Peptides.load_dunbrack(Float64)

    include("Peptides/mutators.jl")
end

