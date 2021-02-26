module Caterpillar

using ProtoSyn
using ProtoSyn.Calculators: EnergyFunctionComponent, full_distance_matrix
using CUDA

"""
    calc_solvation_energy_kernel(coords::CuDeviceArray{T}, Ωi::CuDeviceMatrix{T}, n::Int, rmax::T, sc::T) where {T <: AbstractFloat}

Kernel for CUDA_2 acceleration of `calc_solvation_energy` function.
"""
function calc_solvation_energy_kernel(coords::CuDeviceArray{T}, Ωi::CuDeviceMatrix{T}, n::Int, rmax::T, sc::T) where {T <: AbstractFloat}
    # Note: coords must be in AoS format

    i = ((blockIdx().y - 1) * blockDim().y) + threadIdx().y
    j = ((blockIdx().x - 1) * blockDim().x) + threadIdx().x
    
    if i <= n && j <= n

        _i = i << 2 - (2 + i)
        _j = j << 2 - (2 + j)

        dijx    = coords[_i] - coords[_j]
        dijx_sq = dijx * dijx
        dijy    = coords[_i + 1] - coords[_j + 1]
        dijy_sq = dijy * dijy
        dijz    = coords[_i + 2] - coords[_j + 2]
        dijz_sq = dijz * dijz

        dij_sq = dijx_sq + dijy_sq + dijz_sq
        rij = CUDA.sqrt(dij_sq)

        Ωi[i, j] = 1 - (1 / (1 + exp(sc * (rmax - rij))))
    end

    return
end


"""
    Calculators.calc_solvation_energy([::A], pose::Pose, update_forces::Bool = false; Ω::Int = 24, rmax::T = 12.0, sc::T = 5.0) where {A, T <: AbstractFloat}
    
Calculate the pose solvation energy according to the Caterpillar model. The
model can be fine-tuned: `Ω` defines the minimum number of Cɑ-Cɑ contacts to
consider a residue 'burried'; `rmax` defines the minimum distance between
Cɑ's for a contact to be considered (in Angstrom Å); `sc` defines the 'slope
control' (a higher value defines more sharply when to consider a contact).
Finally, an hydrophobic map can be provided (`hydrophob_map`), with the
penalty/reward for correct burial/exposal of hydrophobic/hydrophilic residues,
respectively. Other map examples can be found in `Peptides.constants.jl`. The
optional `A` parameter defines the acceleration mode used (SISD_0, SIMD_1 or
CUDA_2). If left undefined the default ProtoSyn.acceleration.active mode will be
used. This function does not calculate forces (not applicable), and therefore
the `update_forces` flag serves solely for uniformization with other
energy-calculating functions.

# Examples
```jldoctest
julia> Calculators.calc_solvation_energy(pose)
```
"""
function calc_solvation_energy(::Type{ProtoSyn.CUDA_2}, pose::Pose, update_forces::Bool = false; Ω::Int = 24, rmax::T = 24.0, sc::T = 5.0, hydrophob_map::Dict{String, T} = ProtoSyn.Peptides.doolitle_hydrophobicity_mod7) where {T <: AbstractFloat}
    # coords must be in AoS format
    
    s = (an"CA")(pose, gather = true)
    _size = length(s)
    if _size == 0
        @warn "No CA atoms found in the pose: current solvation energy model is unusable."
        return 0.0, nothing
    end
    coords = zeros(eltype(pose.state), 3, _size)
    for (i, atom) in enumerate(s)
        coords[:, i] = pose.state[atom].t
    end

    _coords = CuArray(coords)
    
    # Define the configuration
    n_threads     = min(_size, 32)
    threads       = (n_threads, n_threads)
    n_blocks      = ceil(Int, _size / n_threads)
    blocks        = (n_blocks, n_blocks)
    
    results = CuArray(zeros(T, _size, _size))
    
    @cuda blocks = blocks threads = threads calc_solvation_energy_kernel(_coords, results, _size, rmax, sc)

    residues = collect(eachresidue(pose.graph))
    e_sol = T(0)

    for i in 1:size(results)[1]
        Ωi = sum(results[:, i])

        dhi = hydrophob_map[residues[i].name]

        if Ωi > Ω
            esol_i = T(0)
        else
            esol_i = dhi * (Ω - Ωi)
        end

        e_sol += esol_i
    end
    
    return e_sol, nothing
end


function calc_solvation_energy(A::Union{Type{ProtoSyn.SISD_0}, Type{ProtoSyn.SIMD_1}}, pose::Pose, update_forces::Bool = false; Ω::Int = 24, rmax::T = 12.0, sc::T = 5.0, hydrophob_map::Dict{String, T} = ProtoSyn.Peptides.doolitle_hydrophobicity_mod7) where {T <: AbstractFloat}
        
    dm       = full_distance_matrix(A, pose, an"CA")
    residues = collect(eachresidue(pose.graph))
    esol = T(0)

    for i in 1:size(dm)[1]
        Ωi = 0
        for j in 1:size(dm)[2]
            Ωi += 1 - (1 / (1 + exp(sc * (rmax - dm[i, j]))))
        end
        dhi = hydrophob_map[residues[i].name]

        if Ωi > Ω
            esol_i = T(0)
        else
            esol_i = dhi * (Ω - Ωi)
        end

        esol += esol_i
    end

    return esol, nothing
end

calc_solvation_energy(pose::Pose, update_forces::Bool = false; Ω::Int = 24, rmax::T = 12.0, sc::T = 5.0, hydrophob_map::Dict{String, T} = ProtoSyn.Peptides.doolitle_hydrophobicity_mod7) where {T <: AbstractFloat} = begin
    return calc_solvation_energy(ProtoSyn.acceleration.active, pose, update_forces; Ω = Ω, rmax = rmax, sc = sc, hydrophob_map = hydrophob_map)
end

# * Default Energy Components --------------------------------------------------

function get_default_solvation_energy(α::T = 1.0) where {T <: AbstractFloat}
    return EnergyFunctionComponent(
        "Caterpillar_Solvation",
        calc_solvation_energy,
        Dict{Symbol, Any}(:Ω => 24, :rmax => 12.0, :sc => 5.0, :hydrophob_map => ProtoSyn.Peptides.doolitle_hydrophobicity_mod7),
        α,
        false)
end
end