module SASA

    using ProtoSyn
    using LinearAlgebra

    """
    # TODO
    """
    function calc_sasa(::Type{A}, pose::Pose, selection::Opt{AbstractSelection} = an"CA", update_forces::Bool = false; probe_radius::T = 1.4, n_points::Int = 100) where {A <: ProtoSyn.AbstractAccelerationType, T <: AbstractFloat}
        
        if selection !== nothing
            sele = ProtoSyn.promote(selection, Atom)
        else
            sele = TrueSelection{Atom}()
        end

        dm = ProtoSyn.Calculators.full_distance_matrix(A, pose, sele)
        if A === ProtoSyn.CUDA_2
            dm = collect(dm)
        end
        
        atoms          = sele(pose, gather = true)
        n_atoms        = length(atoms)
        probe_radius_2 = 2 * probe_radius
        
        # Define sphere
        sphere = ProtoSyn.fibonacci_sphere(T, n_points)

        sasas = Vector{T}()
        clouds = Vector{Vector{T}}()
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
            clouds = vcat(clouds, cloud_i)
            esol += Ωi
            push!(sasas, Ωi)
        end

        return esol, nothing, sasas, clouds
    end

    calc_sasa(pose::Pose, selection::Opt{AbstractSelection} = nothing, update_forces::Bool = false; probe_radius::T = 1.4, n_points::Int = 100) where {T <: AbstractFloat} = begin
        calc_sasa(ProtoSyn.acceleration.active, pose, selection, update_forces, probe_radius = probe_radius, n_points = n_points)
    end


    """
    # TODO
    """
    function get_default_sasa(;α::T = 1.0) where {T <: AbstractFloat}
        return ProtoSyn.Calculators.EnergyFunctionComponent(
            "SASA",
            calc_sasa,
            nothing,
            Dict{Symbol, Any}(
                :n_points           => 100,
                :probe_radius       => 1.4,
            ),
            α,
            false)
    end

    # --------------------------------------------------------------------------

    """
    # TODO
    """
    function calc_sasa_energy(::Type{A}, pose::Pose,
        selection::Opt{AbstractSelection} = an"CA",
        update_forces::Bool = false;
        probe_radius::T = 1.4,
        n_points::Int = 100,
        residue_selection::AbstractSelection = TrueSelection{Residue}(),
        hydrophobicity_map::Opt{Dict{String, T}} = nothing,
        max_sasas::Opt{Dict{String, T}} = nothing,
        Ω::T = 0.0) where {A <: ProtoSyn.AbstractAccelerationType, T <: AbstractFloat}

        @assert 0.0 < Ω < 1.0 "Average exposure value (Ω) must be between 0.0 and 1.0!"
        
        if hydrophobicity_map === nothing; hydrophobicity_map = Dict{String, T}(); end
        if max_sasas === nothing; max_sasas = Dict{String, T}(); end

        esol = T(0.0)

        if selection !== nothing
            sele = ProtoSyn.promote(selection, Atom)
        else
            sele = TrueSelection{Atom}()
        end

        a_resi_ids = [a.container.id for a in sele(pose, gather = true)]

        s, _, sasas = calc_sasa(A, pose, sele, update_forces,
            probe_radius = probe_radius, n_points = n_points)

        resi_sele = ProtoSyn.promote(residue_selection, Residue)
        residues = resi_sele(pose, gather = true)
        for residue in residues
            rname = residue.name
            resi_atoms_indexes = findall(==(residue.id), a_resi_ids)
            resi_sasa = sum(sasas[resi_atoms_indexes])
            σi  = hydrophobicity_map[rname]
            if rname in keys(max_sasas)
                ref = max_sasas[rname]
            else
                @warn "No max SASA entry found for Residue type $rname."
                ref = 0.0
            end
            ProtoSyn.verbose.mode && println("Esol in residue $(residue.id)-$(residue.name): $σi * ($resi_sasa - ($ref * $Ω ($(ref * Ω)))) = $(σi * (resi_sasa - (ref * Ω)))")
            esol += σi * (resi_sasa - (ref * Ω))
        end

        return esol, nothing
    end

    calc_sasa_energy(pose::Pose, selection::Opt{AbstractSelection} = nothing, update_forces::Bool = false; probe_radius::T = 1.4, n_points::Int = 100, residue_selection::AbstractSelection = TrueSelection{Residue}(), hydrophobicity_map::Opt{Dict{String, T}} = nothing, max_sasas::Opt{Dict{String, T}} = nothing, Ω::T = 0.0) where {T <: AbstractFloat} = begin
        calc_sasa_energy(ProtoSyn.acceleration.active, pose, selection, update_forces, probe_radius = probe_radius, n_points = n_points, residue_selection = residue_selection, hydrophobicity_map = hydrophobicity_map, max_sasas = max_sasas, Ω = Ω)
    end
    

    """
    # TODO
    """
    function get_default_sasa_energy(;α::T = 1.0) where {T <: AbstractFloat}
        return ProtoSyn.Calculators.EnergyFunctionComponent(
            "SASA_Solvation",
            calc_sasa_energy,
            nothing,
            Dict{Symbol, Any}(
                :probe_radius                 => 1.4,
                :n_points                     => 100,
                :hydrophobicity_map           => ProtoSyn.Peptides.doolitle_hydrophobicity,
                :reference_energies           => Dict{String, Float64}()
            ),
            α,
            false)
    end

    # --------------------------------------------------------------------------

    """
    # TODO
    """
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


    """
    # TODO
    """
    function calc_radius_gyration_energy(pose::Pose, selection::Opt{AbstractSelection}, update_forces::Bool)
        rg = calc_radius_gyration(pose, selection)
        return sum(rg), nothing
    end


    """
    # TODO
    """
    function get_default_rg(;α::T = 1.0) where {T <: AbstractFloat}
        return ProtoSyn.Calculators.EnergyFunctionComponent(
            "Radius_Gyration",
            calc_radius_gyration_energy,
            nothing,
            Dict{Symbol, Any}(),
            α,
            false)
    end
end