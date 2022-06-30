module SASA

    using ProtoSyn
    using LinearAlgebra

    """
        calc_sasa([::Type{A}], pose::Pose, [selection::Opt{AbstractSelection} = an"CA"], [update_forces::Bool = false]; [probe_radius::T = 1.4], [n_points::Int = 100]) where {A <: ProtoSyn.AbstractAccelerationType, T <: AbstractFloat}

    Approximates the given [`Pose`](@ref) `pose` Surface Accessible Surface Area
    (SASA) using the Overlapping Spheres (OLS) algorithm, a variant of the
    Shrake and Rupley algorithm (for more details, see
    https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2712621/). In this algorithm,
    a sphere of points is generated around each [`Atom`](@ref) instance (using
    the [`fibonacci_sphere`](@ref ProtoSyn.fibonacci_sphere) method). The
    spherical distance to the [`Atom`](@ref) is given by `probe_radius` (1.4 Å,
    by default) and the number of points is given by `n_points` (100, by
    default). For all other considered [`Atom`](@ref) instances, points from the
    generated sphere are removed if they are within `probe_radius` cut-off. The
    SASA value is approximated by the number of remaining points in the
    generated sphere (not within `probe_radius` of any other [`Atom`](@ref)
    instance). If an `AbstractSelection` `selection` is provided, the subset of
    [`Atom`](@ref) instances considered for the calculation is restricted to the
    selected atoms (selects atoms named "CA", by default; will
    [`promote`](@ref ProtoSyn.promote) all provided `AbstractSelection`
    instances to be of [`Atom`](@ref) type). This Calculator does not calculate
    forces. As such, `update_forces` has no effect and exists only in order to
    standardize calls between Calculators. An optional parameter
    `Type{<: AbstractAccelerationType}` can be provided, stating the
    acceleration type used to calculate this energetic contribution (See
    [ProtoSyn acceleration types](@ref), if not provided defaults to
    `ProtoSyn.acceleration.active`). Besides the [`Pose`](@ref) `pose` energy
    and forces (set to `nothing` on this Calculator), also returns the
    individually calculated SASA values (for each considered [`Atom`](@ref)
    instance) and the cartesian coordinates of all remaining points in each
    considered [`Atom`](@ref) sphere.

    # See also
    [`get_default_sasa`](@ref) [`calc_sasa_energy`](@ref)

    # Examples
    ```jldoctest
    julia> ProtoSyn.Calculators.SASA.calc_sasa(pose)
    (36597.0, nothing, [...], [...])
    ```
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
        get_default_sasa(;[α::T = 1.0]) where {T <: AbstractFloat}

    Return the default SASA [`EnergyFunctionComponent`](@ref ProtoSyn.Calculators.EnergyFunctionComponent). `α` sets the
    component weight (on an
    [`EnergyFunction`](@ref ProtoSyn.Calculators.EnergyFunction) instance, `1.0`
    by default). This function employs [`calc_sasa`](@ref) as the `:calc`
    function.

    # Settings
    * `n_points::Int` - The number of points to generate in each [`Atom`](@ref) sphere (higher number of points leads to higher accuracy, at the expense of performance);
    * `probe_radius::T` - The distance of each point in a generated sphere to the central [`Atom`](@ref) instance. Any point within `probe_radius` of any other atom is considered buried [`Residue`](@ref) name (where T <: AbstractFloat).

    # See also
    [`calc_sasa`](@ref)

    # Examples
    ```jldoctest
    julia> ProtoSyn.Calculators.SASA.get_default_sasa()
    🞧  Energy Function Component:
    +---------------------------------------------------+
    | Name           | SASA                             |
    | Alpha (α)      | 1.0                              |
    | Update forces  | false                            |
    | Calculator     | calc_sasa                        |
    +---------------------------------------------------+
     |    +----------------------------------------------------------------------------------+
     ├──  ● Settings                      | Value                                            |
     |    +----------------------------------------------------------------------------------+
     |    | probe_radius                  | 1.4                                              |
     |    | n_points                      | 100                                              |
     |    +----------------------------------------------------------------------------------+
     |    
     └──  ○  Selection: nothing
    ```
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
        calc_sasa_energy([::Type{A}], pose::Pose, [selection::Opt{AbstractSelection} = an"CA"], [update_forces::Bool = false]; [probe_radius::T = 1.4], [n_points::Int = 100], [residue_selection::AbstractSelection = TrueSelection{Residue}()], [hydrophobicity_map::Opt{Dict{String, T}} = nothing], [max_sasas::Opt{Dict{String, T}} = nothing], [Ω::T = 0.0]) where {A <: ProtoSyn.AbstractAccelerationType, T <: AbstractFloat}
    
    Approximates the given [`Pose`](@ref) `pose` Surface Accessible Surface Area
    (SASA) energy using the Overlapping Spheres (OLS) algorithm (See
    [`calc_sasa`](@ref) for more information on how SASA is calculated). In this
    Calculator, the calculated SASA of a [`Residue`](@ref) (using the `n_points`
    sphere at `probe_radius` from all considered [`Atom`](@ref) instances) is
    compared to reference values (`max_sasas`, define the maximum SASA value of
    a [`Residue`](@ref) completly solvated). The average exposure value `Ω`
    defines the average ratio of `max_sasas` that is exposed in a given system.
    Any SASA value bellow the Ω * max_sasa is considered buried, and vice-versa.
    This burial/exposure scale is then multiplied by that [`Residue`](@ref)
    specific hydrophobiciy index (provided in the `hydrophobicity_map`).
    Hydrophilic/exposed and hydrophobic/buried [`Residue`](@ref) instances
    should be rewarded and vice-versa. Since SASA approximation is
    [`Atom`](@ref) based, any provided `AbstractSelection` `selection` limits
    the subset of considered [`Atom`](@ref) instances for SASA calculation. Each
    [`Atom`](@ref) individual contribution is summed to the [`Residue`](@ref)
    level, only for [`Residue`](@ref) instances selected by `residue_selection`.
    This Calculator does not calculate forces. As such, `update_forces` has no
    effect and exists only in order to standardize calls between Calculators. An
    optional parameter `Type{<: AbstractAccelerationType}` can be provided,
    stating the acceleration type used to calculate this energetic contribution
    (See [ProtoSyn acceleration types](@ref), if not provided defaults to
    `ProtoSyn.acceleration.active`).

    # See also
    [`get_default_sasa_energy`](@ref) [`calc_sasa`](@ref)

    # Examples
    ```jldoctest
    julia> ProtoSyn.Calculators.SASA.calc_sasa_energy(pose, hydrophobicity_map = ProtoSyn.Peptides.doolitle_hydrophobicity)
    (-37074.2, nothing)
    ```
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

        @assert 0.0 <= Ω <= 1.0 "Average exposure value (Ω) must be between 0.0 and 1.0!"
        
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
            @info "Esol in residue $(residue.id)-$(residue.name): $σi * ($resi_sasa - ($ref * $Ω ($(ref * Ω)))) = $(σi * (resi_sasa - (ref * Ω)))"
            esol += σi * (resi_sasa - (ref * Ω))
        end

        return esol, nothing
    end

    calc_sasa_energy(pose::Pose, selection::Opt{AbstractSelection} = nothing, update_forces::Bool = false; probe_radius::T = 1.4, n_points::Int = 100, residue_selection::AbstractSelection = TrueSelection{Residue}(), hydrophobicity_map::Opt{Dict{String, T}} = nothing, max_sasas::Opt{Dict{String, T}} = nothing, Ω::T = 0.0) where {T <: AbstractFloat} = begin
        calc_sasa_energy(ProtoSyn.acceleration.active, pose, selection, update_forces, probe_radius = probe_radius, n_points = n_points, residue_selection = residue_selection, hydrophobicity_map = hydrophobicity_map, max_sasas = max_sasas, Ω = Ω)
    end
    

    """
        get_default_sasa_energy(;[α::T = 1.0]) where {T <: AbstractFloat}

    Return the default SASA energy [`EnergyFunctionComponent`](@ref ProtoSyn.Calculators.EnergyFunctionComponent). `α` sets
    the component weight (on an
    [`EnergyFunction`](@ref ProtoSyn.Calculators.EnergyFunction) instance, `1.0`
    by default). This function employs [`calc_sasa`](@ref) as the `:calc`
    function.

    # Settings
    * `n_points::Int` - The number of points to generate in each [`Atom`](@ref) sphere (higher number of points leads to higher accuracy, at the expense of performance);
    * `probe_radius::T` - The distance of each point in a generated sphere to the central [`Atom`](@ref) instance. Any point within `probe_radius` of any other atom is considered buried [`Residue`](@ref) name (where T <: AbstractFloat);
    * `hydrophobicity_map::Dict{String, T}` - A dictionary of hydrophobicity values for each [`Residue`](@ref) name, positive values indicate hydrophobicity and vice-versa (where T <: AbstractFloat);
    * `max_sasas::Dict{String, T}` - A dictionary of max_sasa values (SASA values for fully-solvated [`Residue`](@ref) instances) for each [`Residue`](@ref) name (where T <: AbstractFloat);
    * `Ω::T` - The average exposure value (between 0.0 and 1.0), any SASA value bellow this percentage of max_sasa is considered buried (where T <: AbstractFloat).

    # See also
    [`calc_sasa_energy`](@ref)

    # Examples
    ```jldoctest
    julia> ProtoSyn.Calculators.SASA.get_default_sasa_energy()
    🞧  Energy Function Component:
    +---------------------------------------------------+
    | Name           | SASA_Solvation                   |
    | Alpha (α)      | 1.0                              |
    | Update forces  | false                            |
    | Calculator     | calc_sasa_energy                 |
    +---------------------------------------------------+
     |    +----------------------------------------------------------------------------------+
     ├──  ● Settings                      | Value                                            |
     |    +----------------------------------------------------------------------------------+
     |    | max_sasas                     | Dict{String, Float64}(0 components)              |
     |    | hydrophobicity_map            | Dict{String, Float64}(22 components)             |
     |    | Ω                             | 0.0                                              |
     |    | probe_radius                  | 1.4                                              |
     |    | n_points                      | 100                                              |
     |    +----------------------------------------------------------------------------------+
     |    
     └──  ○  Selection: nothing
    ```
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
                :max_sasas                    => Dict{String, Float64}(),
                :Ω                            => 0.0
            ),
            α,
            false)
    end
end