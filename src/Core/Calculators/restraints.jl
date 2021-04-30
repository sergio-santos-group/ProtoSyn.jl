@doc """
    module Restraints

Contains restraint energy components, such as `bond_distance_restraint`.
"""
module Restraints

    using ProtoSyn
    using ProtoSyn.Calculators: EnergyFunctionComponent

    """
        Calculators.calc_bond_distance_restraint([::Type{A}], pose::Pose, update_forces::Bool = false; x0::T = 2.0) where {A <: ProtoSyn.AbstractAccelerationType, T <: AbstractFloat}
        
    Calculate the [`Pose`](@ref) `pose` bond distance restraint energy according
    to a quadratic potential, based on the cartesian coordinates (make sure the
    [`Pose`](@ref) `pose` is synched, see [`sync!`](@ref)). This potential
    defines a maximum distance for the bond, based on the elements involved in
    the bond (See `ProtoSyn.max_bond_lengths`). If the pair of elements is not
    found, a default `x0` value is used instead. This function iterates over all
    [`Atom`](@ref) instances in the provided [`Pose`](@ref) `pose` (See
    [Counters and Iterators](@ref)) and checks all bonds in each [`Atom`](@ref).
    This means that bonds energy and forces are usually checked twice, as both
    [`Atom`](@ref) instances involved in a bond have records on one another.
    This is usually not a problem, as the energy value is compared between
    frames in a simulation environment and therefore the scale of the value is
    not important. If the `update_forces` flag is set to `true` (`false`, by
    default), also return the calculated forces based on this potential. Note
    that this function assumes [`Atom`](@ref)`.id` entries are synched between
    the [Graph](@ref graph-types) and [State](@ref state-types) (See
    [Indexation](@ref graph-methods-indexation)). An optional parameter
    `Type{<: AbstractAccelerationType}` can be provided, stating the
    acceleration type used to calculate this energetic contribution (See
    [ProtoSyn acceleration types](@ref)).

    !!! ukw "Note:"
        As of ProtoSyn 1.0, this function's acceleration type defaults to
        `SIMD_0` regardless of the requested acceleration type. This may be
        changed in future iterations.

    # See also
    [`get_default_bond_distance_restraint`](@ref)

    # Examples
    ```jldoctest
    julia> Calculators.Restraints.calc_bond_distance_restraint(pose)
    (7.216571959087505, nothing)
    ```
    """
    function calc_bond_distance_restraint(::Type{ProtoSyn.SISD_0}, pose::Pose, update_forces::Bool = false; x0::T = 2.0) where {T <: AbstractFloat}

        e  = 0.0
        if update_forces
            f = zeros(size(pose.state.f))
        end

        for atom in eachatom(pose.graph)
            for bond in atom.bonds
                d = ProtoSyn.distance(pose.state[atom], pose.state[bond])
                
                p = atom.symbol*bond.symbol
                if p in keys(ProtoSyn.Units.max_bond_lengths)
                    _x0 = ProtoSyn.Units.max_bond_lengths[p]
                else
                    println("Pair $p not found in max lengths list.")
                    _x0 = x0
                end

                if d > x0
                    e += (d - _x0) * (d - _x0)
                    if update_forces
                        v = collect(pose.state[atom].t .- pose.state[bond].t)
                        factor1 = 2 * (d - _x0)

                        # * Assumes atom IDs are synched with pose.state
                        f[:, atom.id] -= v * factor1
                        f[:, bond.id] += v * factor1
                    end
                end
            end
        end

        if update_forces
            return e, f
        else
            return e, nothing
        end
    end

    calc_bond_distance_restraint(pose::Pose, update_forces::Bool = false; x0::T = 2.0) where {T <: AbstractFloat} = begin
        calc_bond_distance_restraint(ProtoSyn.acceleration.active, pose, update_forces, x0 = x0)
    end

    calc_bond_distance_restraint(::Type{A}, pose::Pose, update_forces::Bool = false; x0::T = 2.0) where {A <: ProtoSyn.AbstractAccelerationType, T <: AbstractFloat} = begin
        # ! Regardless of the requested acceleration type, this function will
        # ! default to SISD_0. This may be changed in future iterations. 
        calc_bond_distance_restraint(ProtoSyn.SISD_0, pose, update_forces, x0 = x0)
    end


    """
        get_default_bond_distance_restraint(;α::T = 1.0) where {T <: AbstractFloat}

    Return the default bond distance restraint
    [`EnergyFunctionComponent`](@ref). `α` sets the component weight (on an
    [`EnergyFunction`](@ref ProtoSyn.Calculators.EnergyFunction) instance, `1.0` by default). This function employs
    [`calc_bond_distance_restraint`](@ref) as the `:calc` function.

    # Settings
    * `x0::Float64` - The maximum allowed bond distance. Any bond with a longer distance will be subjected to a quadratic energy penalty. This value is normally extracted from `ProtoSyn.Units.max_bond_lengths`. If the pair of [`Atom`](@ref) instances identified in a bond is not found in this table, use this default `x0` value (Default: `2.0`).

    # See also
    [`calc_bond_distance_restraint`](@ref)

    # Examples
    ```jldoctest
    julia> ProtoSyn.Calculators.Restraints.get_default_bond_distance_restraint()
             Name : Bond_Distance_Restraint
       Weight (α) : 1.0
    Update forces : true
          Setings :
             :x0 => 2.0
    ```
    """
    function get_default_bond_distance_restraint(;α::T = 1.0) where {T <: AbstractFloat}
        return EnergyFunctionComponent(
            "Bond_Distance_Restraint",
            calc_bond_distance_restraint,
            Dict{Symbol, Any}(:x0 => 2.0),
            α,
            true)
    end


    # --------------------------------------------------------------------------
    # * Flat bottom restraint base function

    MaskMap = Opt{Union{ProtoSyn.Mask{<: ProtoSyn.AbstractContainer}, Matrix{<: AbstractFloat}, Function}}

    """
        calc_flat_bottom_restraint([::Type{A}], pose::Pose, update_forces::Bool; d1::T = 0.0, d2::T = 0.0, d3::T = Inf, d4::T = Inf, selection::Opt{AbstractSelection} = nothing, mask::MaskMap = nothing) where {A <: ProtoSyn.AbstractAccelerationType, T <: AbstractFloat}

    Apply a flat bottom potential to a given [`Pose`](@ref) `pose`. The
    potential is iteratively obtained each call using the
    [`get_flat_bottom_potential`](@ref ProtoSyn.Calculators.get_flat_bottom_potential)
    method (See [Available potentials](@ref)), by providing the given `d1::T`
    (default: 0.0), `d2::T` (default = 0.0), `d3::T` (default = Inf) and `d4::T`
    (default = Inf) settings as the flat bottom potential distances. This
    potential then applied to the [`Pose`](@ref) `pose` (via the
    [`apply_potential`](@ref ProtoSyn.Calculators.apply_potential) method),
    optionally on a subset of [`Atom`](@ref) instances given by the
    `AbstractSelection` `selection` and optionally multiplied by a `mask`. This
    `mask` can be a [`Mask`](@ref), a `Matrix{T}` or a `Function`, in which case
    it should be a functor (return a `Function`) (For the correct signature of
    this `Function` `mask`, see [Creating custom masks](@ref)). These 3 options
    are named `MaskMap` for a simplicity of organization only. Return the total
    energy of the system and matrix of forces felt on each atom. Note that the
    calculation acceleration type can be set by providing an option parameter
    `Type{<: ProtoSyn.AbstractAccelerationType}`. If not provided, the default
    `ProtoSyn.acceleration.active` will be used instead.

    !!! ukw "Note:"
        As of ProtoSyn 1.0, the
        [`apply_potential`](@ref ProtoSyn.Calculators.apply_potential)
        acceleration type defaults to `CUDA_2` regardless of the requested
        acceleration type. This may be changed in future iterations.

    # Examples
    ```jldoctest
    julia> ProtoSyn.Calculators.Restraints.calc_flat_bottom_restraint(pose, true)
    (0.0, [0.0 0.0 … 0.0 0.0; 0.0 0.0 … 0.0 0.0; 0.0 0.0 … 0.0 0.0])

    julia> ProtoSyn.Calculators.Restraints.calc_flat_bottom_restraint(pose, false, d1 = 10.0, d2 = 12.0)
    41966.43234537284
    ```
    """
    function calc_flat_bottom_restraint(::Type{A}, pose::Pose, update_forces::Bool; d1::T = 0.0, d2::T = 0.0, d3::T = Inf, d4::T = Inf, selection::Opt{AbstractSelection} = nothing, mask::MaskMap = nothing) where {A <: ProtoSyn.AbstractAccelerationType, T <: AbstractFloat}
        fbr = ProtoSyn.Calculators.get_flat_bottom_potential(d1 = d1, d2 = d2, d3 = d3, d4 = d4)
        e, f = ProtoSyn.Calculators.apply_potential(A, pose, fbr, mask, selection)
        if update_forces
            return e, f
        else
            return e
        end
    end # function

    calc_flat_bottom_restraint(pose::Pose, update_forces::Bool; d1::T = 0.0, d2::T = 0.0, d3::T = Inf, d4::T = Inf, selection::Opt{AbstractSelection} = nothing, mask::MaskMap = nothing) where {T <: AbstractFloat} = begin
        calc_flat_bottom_restraint(ProtoSyn.acceleration.active, pose, update_forces, d1 = d1, d2 = d2, d3 = d3, d4 = d4, selection = selection, mask = mask)
    end
end