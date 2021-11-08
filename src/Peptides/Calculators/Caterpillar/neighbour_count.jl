function neighbour_count(::Type{A}, pose::Pose, update_forces::Bool; selection::AbstractSelection = ProtoSyn.TrueSelection{Atom}(), identification_curve::Function = () -> nothing, hydrophobicity_weight::Function = () -> nothing, rmax::T = 9.0, sc::T = 1.0, Ω::Union{Int, T} = 750.0, hydrophobicity_map::Dict{String, T} = ProtoSyn.Peptides.doolitle_hydrophobicity, kwargs...) where {A <: ProtoSyn.AbstractAccelerationType, T <: AbstractFloat}
    
    dm    = ProtoSyn.Calculators.full_distance_matrix(A, pose, selection)
    if A === ProtoSyn.CUDA_2
        dm = collect(dm)
    end
    atoms = selection(pose, gather = true)
    
    # Ωis = Vector{T}()
    # esols = Vector{T}()
    esol = T(0.0)
    for i in 1:size(dm)[1]
        Ωi = 0.0
        for j in 1:size(dm)[2]
            i === j && continue
            Ωi += identification_curve(dm[i, j]; rmax = rmax, sc = sc)
        end

        DHI = hydrophobicity_map[atoms[i].container.name]
        esol_i = hydrophobicity_weight(Ωi, hydrophobicity_map_value = DHI, Ω = Ω)
        esol += esol_i
        
        # push!(Ωis, Ωi)
        # push!(esols, esol_i)
    end

    # return Ωis, esols 
    return esol, nothing
end

neighbour_count(pose::Pose, update_forces::Bool; selection::AbstractSelection = ProtoSyn.TrueSelection{Atom}(), identification_curve::Function = () -> nothing, hydrophobicity_weight::Function = () -> nothing, rmax::T = 9.0, sc::T = 1.0, Ω::Union{Int, T} = 750.0, hydrophobicity_map::Dict{String, T} = ProtoSyn.Peptides.doolitle_hydrophobicity, kwargs...) where {A <: ProtoSyn.AbstractAccelerationType, T <: AbstractFloat} = begin
    neighbour_count(ProtoSyn.acceleration.active, pose, update_forces;
        selection = selection,
        identification_curve = identification_curve,
        hydrophobicity_weight = hydrophobicity_weight,
        rmax = rmax,
        sc = sc, Ω = Ω,
        hydrophobicity_map = hydrophobicity_map)
end