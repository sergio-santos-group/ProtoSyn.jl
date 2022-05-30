
module ProtoSyn
__precompile__()

const _force_recompilation = 1
const _version             = 1.01

include("Core/Methods/log.jl")

if !("JULIA_PROTOSYN_WARN_NON_AVALIABLE_EFC" in keys(ENV))
    ENV["JULIA_PROTOSYN_WARN_NON_AVALIABLE_EFC"] = true
end

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

Base.show(io::IO, accel::Acceleration) = begin
    check = accel.active === ProtoSyn.SISD_0 ? "✓" : ""
    println(io, "SISD_0 $check")
    check = accel.active === ProtoSyn.SIMD_1 ? "✓" : ""
    println(io, "SIMD_1 $check")
    check = accel.active === ProtoSyn.CUDA_2 ? "✓" : ""
    println(io, "CUDA_2 $check")
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

try
    if CUDA.has_cuda() && CUDA.has_cuda_gpu()
        acceleration.active = CUDA_2
    else
        println("CUDA not available on this machine.")
    end
catch LoadError
    @warn "CUDA package not loaded."
end

@info "Current acceleration set to $(acceleration.active)"

const resource_dir = joinpath(dirname(@__DIR__), "resources")

@info "Loading Core"

export Opt
const Opt = Union{Nothing, T} where T

include("Core/Units/Units.jl")
include("Core/Methods/python-packages.jl")
include("Core/Methods/other.jl") # Requires Units.jl
include("Core/Types/constants.jl")
include("Core/Types/other.jl")
include("Core/Methods/macros.jl") # Requires other.jl
include("Core/Types/graph.jl")
include("Core/Types/state.jl")
include("Core/Types/pose.jl")
include("Core/Submodules/Selections/selections.jl")
include("Core/Methods/graph.jl") # Requires selections.jl
include("Core/Methods/measure.jl") # Requires selections.jl
include("Core/Methods/state.jl")
include("Core/Methods/base.jl")
include("Core/Methods/io.jl")
include("Core/Methods/iterators.jl")
include("Core/Methods/align.jl")
include("Core/Submodules/Clustering/Clustering.jl")
include("Core/Submodules/Builder/grammar.jl")
include("Core/Methods/tautomers.jl") # Requires grammar.jl
include("Core/Submodules/Builder/Builder.jl")
include("Core/Methods/pose.jl") # Requires grammar.jl & graph.jl
include("Core/Submodules/Builder/modifications.jl") # Requires pose.jl
include("Core/Submodules/ExternalPackages/GROMACS/gromacs.jl") # Requires pose.jl

@info "Loading Calculators"
include("Core/Calculators/Calculators.jl")
include("Core/Methods/histogram.jl") # Requires distance_matrix.jl

@info "Loading Mutators"
include("Core/Mutators/Mutators.jl")

@info "Loading Drivers"
include("Core/Drivers/Drivers.jl")

@info "Loading Peptides"
include("Peptides/Peptides.jl")

@info "Loading Materials"
include("Materials/Materials.jl")

@info "Loading Sugars"
include("Sugars/Sugars.jl")

@info "Loading Common"
include("Common/Common.jl")

@info "Loading external models & packages" 

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

function __init__()
    @info "ProtoSyn loaded successfully!"
    version()
end

end # module