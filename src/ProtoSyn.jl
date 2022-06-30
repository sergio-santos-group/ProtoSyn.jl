
module ProtoSyn
__precompile__()

const _force_recompilation = 1
const _version             = 1.1

include("Core/Methods/log.jl")

print_loading("External packages")
using CpuId
using Printf

printstyled(" | Loading SIMD\n", color = :cyan)
using SIMD
printstyled(" | Loading CUDA\n", color = :cyan)
using CUDA

print_loading("Setting up global variables")
abstract type AbstractAccelerationType end
abstract type SISD_0 <: AbstractAccelerationType end
abstract type SIMD_1 <: AbstractAccelerationType end
abstract type CUDA_2 <: AbstractAccelerationType end

export Opt
const Opt = Union{Nothing, T} where T

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

printstyled(" | Current acceleration set to $(acceleration.active)\n", color = :cyan)

const resource_dir = joinpath(dirname(@__DIR__), "resources")

print_loading("Core module")
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

printstyled(" | Loading Calculators\n", color = :cyan)
include("Core/Calculators/Calculators.jl")
include("Core/Methods/histogram.jl") # Requires distance_matrix.jl

printstyled(" | Loading Mutators\n", color = :cyan)
include("Core/Mutators/Mutators.jl")

printstyled(" | Loading Drivers\n", color = :cyan)
include("Core/Drivers/Drivers.jl")

print_loading("Peptides module")
include("Peptides/Peptides.jl")

print_loading("Materials module")
include("Materials/Materials.jl")

print_loading("Sugars module")
include("Sugars/Sugars.jl")

print_loading("Common module")
include("Common/Common.jl")

print_loading("External models")

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
    println("\n")
end

function __init__()
    print_loading("ProtoSyn loaded successfully!", color = :green)
    version()
    set_logger_to_warn()
end

end # module