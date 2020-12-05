module ProtoSyn

# __precompile__(true)

@info "Loading required packages"
using CpuId
@info " | Loading SIMD"
@time using SIMD
@info " | Loading CUDA"
@time using CUDA

@info "Setting up variables"
@time begin

    abstract type AbstractAccelerationType end
    abstract type SISD_0 <: AbstractAccelerationType end
    abstract type SIMD_1 <: AbstractAccelerationType end
    abstract type CUDA_2 <: AbstractAccelerationType end

    mutable struct Acceleration
        active::Type{<: AbstractAccelerationType}
    end

    acceleration = Acceleration(SISD_0)
    try
        if simdbits() >= 256
            acceleration.active = SIMD_1
        else
            println("SIMD not available on this machine.")
        end
    catch LoadError
        @warn "SIMD package not loaded."
    end

    cuda_available = false
    try
        if CUDA.has_cuda() && CUDA.has_cuda_gpu()
            acceleration.active = CUDA_2
        else
            println("CUDA not available on this machine.")
        end
    catch LoadError
        @warn "CUDA package not loaded."
    end

    println("Current acceleration set to $acceleration")
end

const resource_dir = joinpath(dirname(@__DIR__), "resources")

@info "Loading Core"
@time begin
    include("Core/XMLRPC/XMLRPC.jl")
    include("Core/Units/Units.jl")
    include("Core/Methods/constants.jl")
    include("Core/Methods/macros.jl")
    include("Core/Types/graph.jl")
    include("Core/Types/state.jl")
    include("Core/Types/pose.jl")
    include("Core/Methods/graph.jl")
    include("Core/Methods/measure.jl")
    include("Core/Methods/state.jl")
    include("Core/Methods/pose.jl")
    include("Core/Methods/base.jl")
    include("Core/Methods/io.jl")
    include("Core/Methods/iterators.jl")
    include("Core/Selections/selections.jl")
    include("Core/Methods/aux.jl")
    include("Core/Builder/grammar.jl")
    include("Core/Builder/Builder.jl")
end

@info "Loading Calculators"
include("Core/Calculators/Calculators.jl")

@info "Loading Mutators"
include("Core/Mutators/Mutators.jl")

@info "Loading Peptides"
@time include("Peptides/Peptides.jl")

@info "Loading Drivers"
include("Drivers/Drivers.jl")

@info "Loading Common"
include("Common/Common.jl")


@info "ProtoSyn loaded successfully!"
end # module

