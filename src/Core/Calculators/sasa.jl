module SASA

    using ProtoSyn
    using LinearAlgebra

    function calc_exposed_sasa(::Type{A}, pose::Pose, update_forces::Bool = false; selection::AbstractSelection = an"CA", probe_radius::T = 6.0, n_points::Int = 100, hydrophobicity_map::Dict{String, T} = ProtoSyn.Peptides.doolitle_hydrophobicity) where {A <: ProtoSyn.AbstractAccelerationType, T <: AbstractFloat}
        dm = ProtoSyn.Calculators.full_distance_matrix(A, pose, selection)
        if A === ProtoSyn.CUDA_2
            dm = collect(dm)
        end
        
        atoms = selection(pose, gather = true)
        probe_radius_2 = 2 * probe_radius
        
        # Define sphere
        sphere = ProtoSyn.fibonacci_sphere(T, n_points)

        # Ωis   = Vector{T}() # !
        # esols = Vector{T}() # !
        esol  = T(0.0)
        for i in 1:size(dm)[1]
            Ωi = 0

            # define cloud i
            cloud_i  = copy(sphere)
            i_coords = pose.state[atoms[i]].t
            cloud_i  = map(x -> (x .* probe_radius) .+ i_coords, cloud_i)

            for j in 1:size(dm)[2]
                if (dm[i, j] > probe_radius_2) | (i === j)
                    continue
                end

                # measure distances between atom j and each point of the cloud i
                j_coords = pose.state[atoms[j]].t
                for k in length(cloud_i):-1:1 # note the reverse loop; can't use n_points -> length(cloud_i) should be changing
                    point_i = cloud_i[k]
                    if norm(point_i - j_coords) < probe_radius
                        deleteat!(cloud_i, k)
                    end
                end
            end

            Ωi = length(cloud_i)
            esol_i = Ωi * hydrophobicity_map[atoms[i].container.name]
            esol += esol_i
            
            # push!(Ωis, Ωi) # !
            # push!(esols, esol_i) # !
        end

        return esol, nothing
        # return esol, nothing, Ωis, esols 
    end

    calc_exposed_sasa(pose::Pose, update_forces::Bool = false; selection::AbstractSelection = an"CA", probe_radius::T = 6.0, n_points::Int = 100, hydrophobicity_map::Dict{String, T} = ProtoSyn.Peptides.doolitle_hydrophobicity) where {T <: AbstractFloat} = begin
        calc_exposed_sasa(ProtoSyn.acceleration.active, pose, update_forces, selection = selection, probe_radius = probe_radius, n_points = n_points, hydrophobicity_map = hydrophobicity_map)
    end

    # ---

    function calc_core_sasa(::Type{A}, pose::Pose, update_forces::Bool = false; selection::AbstractSelection = an"CA", probe_radius::T = 6.0, n_points::Int = 1000, hydrophobicity_map::Dict{String, T} = ProtoSyn.Peptides.doolitle_hydrophobicity) where {A <: ProtoSyn.AbstractAccelerationType, T <: AbstractFloat}
        dm = ProtoSyn.Calculators.full_distance_matrix(A, pose, selection)
        if A === ProtoSyn.CUDA_2
            dm = collect(dm)
        end
        
        atoms = selection(pose, gather = true)
        probe_radius_2 = 2 * probe_radius
        
        # Define sphere
        cloud  = ProtoSyn.fibonacci_sphere(T, n_points)
        levels = repeat([T(0.0)], n_points)

        # _cloud  = Vector{Vector{T}}()
        # _levels = Vector{T}()

        esol = T(0.0)
        
        for i in 1:size(dm)[1]
            # define cloud i
            cloud_i  = copy(cloud)
            levels_i = copy(levels)
            i_coords = pose.state[atoms[i]].t
            cloud_i  = map(x -> (x .* probe_radius) .+ i_coords, cloud_i)

            for j in i:size(dm)[2]
                if (dm[i, j] > probe_radius_2)
                    continue
                end

                j_coords = pose.state[atoms[j]].t
                for k in 1:length(cloud_i) # note the reverse loop; can't use n_points -> length(cloud_i) should be changing
                    if norm(cloud_i[k] - j_coords) < probe_radius
                        levels_i[k] += hydrophobicity_map[atoms[j].container.name]
                    end
                end
            end

            # for point in cloud_i
            #     push!(_cloud, point)
            # end

            # for level in levels_i
            #     push!(_levels, level)
            # end

            # Define what to count
            esol += sum(filter(x -> x < 0.0, levels_i))
        end

        # return esol, nothing, _cloud, _levels
        return esol, nothing
    end

    calc_core_sasa(pose::Pose, update_forces::Bool = false; selection::AbstractSelection = an"CA", probe_radius::T = 6.0, n_points::Int = 1000, hydrophobicity_map::Dict{String, T} = ProtoSyn.Peptides.doolitle_hydrophobicity) where {T <: AbstractFloat} = begin
        calc_core_sasa(ProtoSyn.acceleration.active, pose, update_forces, selection = selection, probe_radius = probe_radius, n_points = n_points, hydrophobicity_map = hydrophobicity_map)
    end

    function get_default_sasa(;α::T = 1.0) where {T <: AbstractFloat}
        return ProtoSyn.Calculators.EnergyFunctionComponent(
            "SASA_CORE_solvation",
            calc_core_sasa,
            Dict{Symbol, Any}(
                :selection          => an"CA",
                :probe_radius       => 6.0,
                :n_points           => 1000,
                :hydrophobicity_map => ProtoSyn.Peptides.doolitle_hydrophobicity
            ),
            α,
            false
        )
    end

end