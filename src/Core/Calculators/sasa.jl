module SASA

    using ProtoSyn
    using LinearAlgebra

    function calc_sasa(::Type{A}, pose::Pose, selection::AbstractSelection = an"CA", update_forces::Bool = false; probe_radius::T = 6.0, n_points::Int = 100, hydrophobicity_map::Dict{String, T} = ProtoSyn.Peptides.doolitle_hydrophobicity, Ω::T = 4.0) where {A <: ProtoSyn.AbstractAccelerationType, T <: AbstractFloat}
        dm = ProtoSyn.Calculators.full_distance_matrix(A, pose, selection)
        if A === ProtoSyn.CUDA_2
            dm = collect(dm)
        end
        
        atoms = selection(pose, gather = true)
        probe_radius_2 = 2 * probe_radius
        
        # Define sphere
        sphere = ProtoSyn.fibonacci_sphere(T, n_points)

        Ωis   = Vector{T}() # !
        esols = Vector{T}() # !
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

            Ωi = length(cloud_i) - Ω
            esol_i = Ωi * hydrophobicity_map[atoms[i].container.name]
            esol += esol_i
            
            push!(Ωis, Ωi) # !
            push!(esols, esol_i) # !
        end

        # return esol, nothing
        return esol, nothing, Ωis, esols 
    end

    calc_sasa(pose::Pose, selection::AbstractSelection = an"CA", update_forces::Bool = false; probe_radius::T = 6.0, n_points::Int = 100, hydrophobicity_map::Dict{String, T} = ProtoSyn.Peptides.doolitle_hydrophobicity, Ω::T = 4.0) where {T <: AbstractFloat} = begin
        calc_sasa(ProtoSyn.acceleration.active, pose, update_forces, selection = selection, probe_radius = probe_radius, n_points = n_points, hydrophobicity_map = hydrophobicity_map, Ω = Ω)
    end

    function get_default_sasa(;α::T = 1.0) where {T <: AbstractFloat}
        return ProtoSyn.Calculators.EnergyFunctionComponent(
            "SASA",
            calc_sasa,
            an"CA",
            Dict{Symbol, Any}(
                :Ω                  => 4.0,
                :n_points           => 100,
                :probe_radius       => 6.0,
                :hydrophobicity_map => ProtoSyn.Peptides.doolitle_hydrophobicity,
            ),
            α,
            false
        )
    end

    # --------------------------------------------------------------------------

    function calc_sasa_solvation_energy(::Type{A}, pose::Pose, update_forces::Bool = false;
        solvation_selection::AbstractSelection = an"CA",
        solvation_probe_radius::T = 5.0, solvation_n_points::Int = 10,
        solvation_cut_off_radius::T = 3.0,
        selection::AbstractSelection = an"CA", probe_radius::T = 6.0,
        clash_radius::T = 3.0, n_points::Int = 1000,
        hydrophobicity_map::Dict{String, T} = ProtoSyn.Peptides.doolitle_hydrophobicity_mod7,
        solvent_hydrophobicity_index::T = -2.0,
        pose_solvated::Opt{Pose} = nothing) where {A <: ProtoSyn.AbstractAccelerationType, T <: AbstractFloat} 
    
        # Calculate new solvent positions
        if pose_solvated === nothing
            pose_solvated = ProtoSyn.add_semi_explicit_solvent(pose,
                selection = solvation_selection,
                probe_radius = solvation_probe_radius,
                n_points = solvation_n_points,
                cut_off_radius = solvation_cut_off_radius)
        end
    
        @assert count(an"SOL"(pose_solvated)) > 0 "No solvent SOL found in the provided pose"
    
        _hydrophobicity_map = copy(hydrophobicity_map)
        _hydrophobicity_map["SOL"] = solvent_hydrophobicity_index
    
        _selection = selection | an"SOL"
        dm = ProtoSyn.Calculators.full_distance_matrix(A, pose_solvated, _selection)
        if A === ProtoSyn.CUDA_2
            dm = collect(dm)
        end
        
        atoms = _selection(pose_solvated, gather = true)
        M = length(atoms)
        N = count(selection(pose_solvated))
        probe_radius_2 = 2 * probe_radius
        
        # Define sphere
        cloud  = ProtoSyn.fibonacci_sphere(T, n_points)
        levels = repeat([T(0.0)], n_points)
    
        _cloud  = Vector{Vector{T}}()
        _levels = Vector{T}()
    
        esol = T(0.0)
        
        for i in 1:N
    
            i_hydro_index = abs(_hydrophobicity_map[atoms[i].container.name])
            i_hydrophobic = _hydrophobicity_map[atoms[i].container.name] >= 0.0
    
            # define cloud i
            cloud_i  = copy(cloud)
            levels_i = copy(levels)
            i_coords = pose_solvated.state[atoms[i]].t
            cloud_i  = map(x -> (x .* probe_radius) .+ i_coords, cloud_i)
    
            for j in 1:M
                if i === j || (dm[i, j] > probe_radius_2)
                    continue
                end
    
                j_hydro_index = abs(_hydrophobicity_map[atoms[j].container.name])
                j_hydrophobic = _hydrophobicity_map[atoms[j].container.name] > 0.0
    
                m = i_hydrophobic === j_hydrophobic ? 1.0 : -1.0
    
                j_coords = pose_solvated.state[atoms[j]].t
                for k in length(cloud_i):-1:1
                    d = norm(cloud_i[k] - j_coords)
                    if (i - 4) < j < (i + 4) && d < clash_radius
                        deleteat!(cloud_i, k)
                        deleteat!(levels_i, k)
                    elseif d < probe_radius
                        levels_i[k] += (i_hydro_index + j_hydro_index) * m
                    end
                end
            end
    
            for point in cloud_i
                push!(_cloud, point)
            end
    
            for level in levels_i
                push!(_levels, level)
            end
    
            # Define what to count
            esol -= sum(levels_i) 
        end
    
        return esol, nothing, _cloud, _levels, pose_solvated
        # return esol, nothing
    end
    
    calc_sasa_solvation_energy(pose::Pose, update_forces::Bool = false;
        solvation_selection::AbstractSelection = an"CA",
        solvation_probe_radius::T = 5.0, solvation_n_points::Int = 10,
        solvation_cut_off_radius::T = 3.0,
        selection::AbstractSelection = an"CA", probe_radius::T = 6.0,
        clash_radius::T = 3.0, n_points::Int = 1000,
        hydrophobicity_map::Dict{String, T} = ProtoSyn.Peptides.doolitle_hydrophobicity_mod7,
        solvent_hydrophobicity_index::T = -2.0,
        pose_solvated::Opt{Pose} = nothing) where {T <: AbstractFloat} = begin

        calc_sasa_solvation_energy(ProtoSyn.acceleration.active, pose, update_forces,
            selection = selection, probe_radius = probe_radius,
            clash_radius = clash_radius, n_points = n_points,
            hydrophobicity_map = hydrophobicity_map,
            solvent_hydrophobicity_index = solvent_hydrophobicity_index,
            solvation_selection = solvation_selection,
            solvation_probe_radius = solvation_probe_radius,
            solvation_n_points = solvation_n_points,
            solvation_cut_off_radius = solvation_cut_off_radius,
            pose_solvated = pose_solvated)
    end

    # ! RE-DO SASA SOLVATION !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    # function get_default_sasa(;α::T = 1.0) where {T <: AbstractFloat}
    #     return ProtoSyn.Calculators.EnergyFunctionComponent(
    #         "SASA_SOL_solvation",
    #         calc_sasa_solvation_energy,
    #         Dict{Symbol, Any}(
    #             :solvation_selection          => an"CA",
    #             :solvation_probe_radius       => 10.0,
    #             :solvation_n_points           => 10,
    #             :solvation_cut_off_radius     => 3.0,
    #             :selection                    => an"CA",
    #             :probe_radius                 => 6.0,
    #             :n_points                     => 100,
    #             :clash_radius                 => 3.0,
    #             :solvent_hydrophobicity_index => -1.0,
    #             :hydrophobicity_map           => ProtoSyn.Peptides.doolitle_hydrophobicity,
    #             :pose_solvated                => nothing
    #         ),
    #         α,
    #         false
    #     )
    # end

    # --------------------------------------------------------------------------

    function calc_radius_gyration(pose::Pose, selection::Opt{AbstractSelection})
        if selection === nothing
            sele = TrueSelection{Atom}()
        else
            sele = ProtoSyn.promote(selection, Atom)
        end

        atoms = sele(pose, gather = true)

        cm = ProtoSyn.center_of_mass(pose, selection)
        mi = [ProtoSyn.Units.mass[atom.symbol] for atom in atoms]'
        coords = pose.state.x.coords[:, sele(pose).content]
        rr = mi .* ((coords .- cm) .^ 2)
        return sqrt.(sum(rr, dims = 2) ./ sum(mi))
    end

    calc_radius_gyration(pose::Pose) = calc_radius_gyration(pose, nothing)

    function calc_radius_gyration_energy(pose::Pose, selection::Opt{AbstractSelection}, update_forces::Bool)
        rg = calc_radius_gyration(pose, selection)
        return sum(rg), nothing
    end

    function get_default_rg(;α::T = 1.0) where {T <: AbstractFloat}
        return ProtoSyn.Calculators.EnergyFunctionComponent(
            "Radius_Gyration",
            calc_radius_gyration_energy,
            selection,
            Dict{Symbol, Any}(),
            α,
            false
        )
    end
end