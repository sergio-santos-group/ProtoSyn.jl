module ProtoSyn

const _version = 0.4

@info "Loading required packages"
using CpuId
using Printf

@info " | Loading SIMD"
using SIMD
@info " | Loading CUDA"
using CUDA

@info "Setting up variables"

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

@info "Current acceleration set to $acceleration"

const resource_dir = joinpath(dirname(@__DIR__), "resources")

@info "Loading Core"

include("Core/XMLRPC/XMLRPC.jl")
include("Core/Units/Units.jl")
include("Core/Methods/constants.jl")
include("Core/Methods/macros.jl")
include("Core/Types/other.jl")
include("Core/Types/graph.jl")
include("Core/Types/state.jl")
include("Core/Types/pose.jl")
include("Core/Submodules/Selections/selections.jl")
include("Core/Methods/graph.jl") # After selections.jl
include("Core/Methods/measure.jl") # After selections.jl
include("Core/Methods/state.jl")
include("Core/Methods/pose.jl")
include("Core/Methods/base.jl")
include("Core/Methods/io.jl")
include("Core/Methods/iterators.jl")
include("Core/Methods/aux.jl")
include("Core/Methods/align.jl")
include("Core/Clustering/Clustering.jl")
include("Core/Submodules/Builder/grammar.jl")
include("Core/Submodules/Builder/Builder.jl")

@info "Loading Calculators"
include("Core/Calculators/Calculators.jl")

abstract type DriverState end
abstract type Driver end

@info "Loading Mutators"
include("Core/Mutators/Mutators.jl")

@info "Loading Drivers"
include("Core/Drivers/Drivers.jl")

@info "Loading Peptides"
include("Peptides/Peptides.jl")

@info "Loading Materials"
include("Materials/Materials.jl")

@info "Loading Common"
include("Common/Common.jl")



@info "ProtoSyn loaded successfully!"

function version()
    header = """
.      ____            _       ____              
      |  _ \\ _ __ ___ | |_ ___/ ___| _   _ _ __  
      | |_) | '__/ _ \\| __/ _ \\___ \\| | | | '_ \\ 
      |  __/| | | (_) | || (_) |__) | |_| | | | |
      |_|   |_|  \\___/ \\__\\___/____/ \\__, |_| |_|
                                       |_/       
    """
    println("\n", header)
    @printf("      %s\n\n", "-"^45)
    @printf(" Version      : %4.2f\n", ProtoSyn._version)
    println(" License      : GNU-GPL-3")
    println(" Developed by : José Pereira (jose.manuel.pereira@ua.pt)")
    println("                Sérgio Santos")
end

version()

end # module

