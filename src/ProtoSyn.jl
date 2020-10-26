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

@info "Loading Peptides"
@time include("Peptides/Peptides.jl")

@info "Loading Peptides"
@time include("Peptides/Peptides.jl")

@info "Loading Calculators"
include("Calculators/Calculators.jl")

@info "ProtoSyn loaded successfully!"
end # module

