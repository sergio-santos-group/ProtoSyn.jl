@doc """
    module Restraints

Contains restraint energy components, such as `bond_distance_restraint`.
"""
module Restraints

    using ProtoSyn
    using ProtoSyn.Calculators: EnergyFunctionComponent, VerletList, MaskMap

    """
        calc_bond_distance_restraint([::Type{A}], pose::Pose, update_forces::Bool = false; x0::T = 2.0) where {A <: ProtoSyn.AbstractAccelerationType, T <: AbstractFloat}
        
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
    [Indexation](@ref core-graph-methods-indexation)). An optional parameter
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
    julia> ProtoSyn.Calculators.Restraints.calc_bond_distance_restraint(pose)
    (0.0, nothing)
    ```
    """
    function calc_bond_distance_restraint(::Type{ProtoSyn.SISD_0}, pose::Pose, selection::Opt{AbstractSelection} = nothing, update_forces::Bool = false; x0::T = 2.0) where {T <: AbstractFloat}

        if selection === nothing
            sele = TrueSelection{Atom}()
        else
            sele = ProtoSyn.promote(selection, Atom)
        end
        
        mask = sele(pose)

        e  = 0.0
        if update_forces
            f = zeros(3, pose.state.size)
        end

        for atom in eachatom(pose.graph)

            !mask[atom.index] && continue

            for bond in atom.bonds

                !mask[bond.index] && continue

                d = ProtoSyn.distance(pose.state[atom], pose.state[bond])
                
                p = atom.symbol*bond.symbol
                if p in keys(ProtoSyn.Units.max_bond_lengths)
                    _x0 = ProtoSyn.Units.max_bond_lengths[p]
                else
                    @info "Pair $p not found in max lengths list."
                    _x0 = x0
                end

                if d > _x0
                    e += (d - _x0) * (d - _x0)
                    
                    if update_forces
                        v = collect(pose.state[atom].t .- pose.state[bond].t)
                        factor1 = 2 * (d - _x0)

                        # * Assumes atom IDs are synched with pose.state
                        f[:, atom.index] -= v * factor1
                        f[:, bond.index] += v * factor1
                    end
                end
            end
        end

        if update_forces
            return e, f[:, mask.content]
        else
            return e, nothing
        end
    end

    calc_bond_distance_restraint(pose::Pose, selection::Opt{AbstractSelection} = nothing, update_forces::Bool = false; x0::T = 2.0) where {T <: AbstractFloat} = begin
        calc_bond_distance_restraint(ProtoSyn.acceleration.active, pose, selection, update_forces, x0 = x0)
    end

    calc_bond_distance_restraint(::Type{A}, pose::Pose, selection::Opt{AbstractSelection} = nothing, update_forces::Bool = false; x0::T = 2.0) where {A <: ProtoSyn.AbstractAccelerationType, T <: AbstractFloat} = begin
        # ! Regardless of the requested acceleration type, this function will
        # ! default to SISD_0. This may be changed in future iterations. 
        calc_bond_distance_restraint(ProtoSyn.SISD_0, pose, selection, update_forces, x0 = x0)
    end


    """
        get_default_bond_distance_restraint(;[α::T = 1.0]) where {T <: AbstractFloat}

    Return the default bond distance restraint
    [`EnergyFunctionComponent`](@ref). `α` sets the component weight (on an
    [`EnergyFunction`](@ref ProtoSyn.Calculators.EnergyFunction) instance, `1.0`
    by default). This function employs
    [`calc_bond_distance_restraint`](@ref)
    as the `:calc` function.

    # Settings
    * `x0::Float64` - The maximum allowed bond distance. Any bond with a longer distance will be subjected to a quadratic energy penalty. This value is normally extracted from `ProtoSyn.Units.max_bond_lengths`. If the pair of [`Atom`](@ref) instances identified in a bond is not found in this table, use this default `x0` value (Default: `2.0`).

    # See also
    [`calc_bond_distance_restraint`](@ref)

    # Examples
    ```jldoctest
    julia> ProtoSyn.Calculators.Restraints.get_default_bond_distance_restraint()
    🞧  Energy Function Component:
    +---------------------------------------------------+
    | Name           | Bond_Distance_Rest               |
    | Alpha (α)      | 1.0                              |
    | Update forces  | true                             |
    | Calculator     | calc_bond_distance_restraint     |
    +---------------------------------------------------+
     |    +----------------------------------------------------------------------------------+
     ├──  ● Settings                      | Value                                            |
     |    +----------------------------------------------------------------------------------+
     |    | x0                            | 2.0                                              |
     |    +----------------------------------------------------------------------------------+
     |    
     └──  ○  Selection: nothing
    ```
    """
    function get_default_bond_distance_restraint(;α::T = 1.0) where {T <: AbstractFloat}
        return EnergyFunctionComponent(
            "Bond_Distance_Rest",
            calc_bond_distance_restraint,
            nothing,
            Dict{Symbol, Any}(:x0 => 2.0),
            α,
            true)
    end

    # --------------------------------------------------------------------------
    # * Flat bottom restraint base function

    """
        calc_flat_bottom_restraint([::Type{A}], pose::Pose, update_forces::Bool; d1::T = 0.0, d2::T = 0.0, d3::T = Inf, d4::T = Inf, selection::Opt{AbstractSelection} = nothing, mask::MaskMap = nothing) where {A <: ProtoSyn.AbstractAccelerationType, T <: AbstractFloat}

    Apply a flat bottom potential to a given [`Pose`](@ref) `pose`. The
    potential is iteratively obtained each call using the
    [`get_flat_bottom_potential`](@ref ProtoSyn.Calculators.get_flat_bottom_potential)
    method (See [Available potentials](@ref)), by providing the given `d1::T`
    (default: 0.0), `d2::T` (default = 0.0), `d3::T` (default = Inf) and `d4::T`
    (default = Inf) settings as the flat bottom potential distances. This
    potential then applied to the [`Pose`](@ref) `pose` (via the
    [`apply_potential!`](@ref ProtoSyn.Calculators.apply_potential!) method),
    optionally on a subset of [`Atom`](@ref) instances given by the
    `AbstractSelection` `selection` and optionally multiplied by a `mask`. This
    `mask` can be a [`Mask`](@ref), a `Matrix{T}` or a `Function`, in which case
    it should be a functor (return a `Function`) (For the correct signature of
    this `Function` `mask`, see [Creating custom masks](@ref)). These 3 options
    are named `MaskMap` for a simplicity of organization only. Return the total
    energy of the system and matrix of forces felt on each atom. Note that the
    calculation acceleration type can be set by providing an option parameter
    `A` `Type{<: ProtoSyn.AbstractAccelerationType}`. If not provided, the
    default `ProtoSyn.acceleration.active` will be used instead.

        calc_flat_bottom_restraint!([::Type{A}], pose::Pose, update_forces::Bool; d1::T = 0.0, d2::T = 0.0, d3::T = Inf, d4::T = Inf, selection::Opt{AbstractSelection} = nothing, mask::MaskMap = nothing) where {A <: ProtoSyn.AbstractAccelerationType, T <: AbstractFloat}

    Apply a flat bottom potential to a given [`Pose`](@ref) `pose` (see above).
    Also apply any energy and forces changes directly to the [`Pose`](@ref)
    `pose`.


    # Examples
    ```
    julia> ProtoSyn.Calculators.Restraints.calc_flat_bottom_restraint(pose, true)
    (0.0, [0.0 0.0 … 0.0 0.0; 0.0 0.0 … 0.0 0.0; 0.0 0.0 … 0.0 0.0])

    julia> ProtoSyn.Calculators.Restraints.calc_flat_bottom_restraint(pose, false, d1 = 10.0, d2 = 12.0)
    (556449.1936070402, [-711.7603616347209 -630.2662235401388 … 995.0284325254745 1153.572133762037; -419.1275359380875 -548.0506257124055 … 286.5285847489888 92.16862928705675; 6.007398880372552 8.2409631821887 … -99.38257889245355 -92.37110004070036])    ```
    ```
    """
    function calc_flat_bottom_restraint(::Type{A}, pose::Pose, selection::Opt{AbstractSelection}, update_forces::Bool; d1::T = 0.0, d2::T = 0.0, d3::T = Inf, d4::T = Inf, mask::MaskMap = nothing, vlist::Opt{VerletList} = nothing) where {A <: ProtoSyn.AbstractAccelerationType, T <: AbstractFloat}
        fbr = ProtoSyn.Calculators.get_flat_bottom_potential(A; d1 = d1, d2 = d2, d3 = d3, d4 = d4)
        e, f = ProtoSyn.Calculators.apply_potential(A, pose, fbr, update_forces, vlist, selection, mask)
        return e, f
    end # function

    calc_flat_bottom_restraint(pose::Pose, selection::Opt{AbstractSelection}, update_forces::Bool; d1::T = 0.0, d2::T = 0.0, d3::T = Inf, d4::T = Inf, mask::MaskMap = nothing, vlist::Opt{VerletList} = nothing) where {T <: AbstractFloat} = begin
        calc_flat_bottom_restraint(ProtoSyn.acceleration.active, pose, selection, update_forces, d1 = d1, d2 = d2, d3 = d3, d4 = d4, mask = mask, vlist = vlist)
    end

    function calc_flat_bottom_restraint!(::Type{A}, pose::Pose, selection::Opt{AbstractSelection}, update_forces::Bool; d1::T = 0.0, d2::T = 0.0, d3::T = Inf, d4::T = Inf, mask::MaskMap = nothing, vlist::Opt{VerletList} = nothing) where {A <: ProtoSyn.AbstractAccelerationType, T <: AbstractFloat}
        fbr = ProtoSyn.Calculators.get_flat_bottom_potential(d1 = d1, d2 = d2, d3 = d3, d4 = d4)
        e, f = ProtoSyn.Calculators.apply_potential!(A, pose, fbr, update_forces, vlist, selection, mask)
    end # function

    calc_flat_bottom_restraint!(pose::Pose, selection::Opt{AbstractSelection}, update_forces::Bool; d1::T = 0.0, d2::T = 0.0, d3::T = Inf, d4::T = Inf, mask::MaskMap = nothing, vlist::Opt{VerletList} = nothing) where {T <: AbstractFloat} = begin
        calc_flat_bottom_restraint!(ProtoSyn.acceleration.active, pose, selection, update_forces, d1 = d1, d2 = d2, d3 = d3, d4 = d4, mask = mask, vlist = vlist)
    end


    """
        get_default_all_atom_clash_restraint(;[α::T = 1.0]) where {T <: AbstractFloat}

    Return the all-atom clash [`EnergyFunctionComponent`](@ref). `α` sets the
    component weight (on an
    [`EnergyFunction`](@ref ProtoSyn.Calculators.EnergyFunction) instance). This
    component employs the [`calc_flat_bottom_restraint`](@ref) method, therefore
    defining a [`Pose`](@ref) energy based on a flat-bottom potential function
    applied to all atom-pairs in the system (N² complexity). By default, this
    [`EnergyFunctionComponent`](@ref) potential sets `:d1` and `:d2` of the
    flat-bottom potential to be 1.0 and 2.0, and masks out bonded atom-pairs.

    # See also
    [`get_default_bond_distance_restraint`](@ref)

    # Settings
    * `d1::T` - The :d1 distance in the flat-bottom potential;
    * `d2::T` - The :d2 distance in the flat-bottom potential;
    * `d3::T` - The :d3 distance in the flat-bottom potential;
    * `d4::T` - The :d4 distance in the flat-bottom potential;
    * `mask::MaskMap` - The [`Mask`](@ref), `Matrix{T}` or `Function` (see [Creating custom masks](@ref)) that masks out our multiplied by a set of pre-defined [`Atom`](@ref) instances;

    # Examples
    ```jldoctest
    julia> ProtoSyn.Calculators.Restraints.get_default_all_atom_clash_restraint()
    🞧  Energy Function Component:
    +---------------------------------------------------+
    | Name           | All_Atom_Clash_Rest              |
    | Alpha (α)      | 1.0                              |
    | Update forces  | true                             |
    | Calculator     | calc_flat_bottom_restraint       |
    +---------------------------------------------------+
     |    +----------------------------------------------------------------------------------+
     ├──  ● Settings                      | Value                                            |
     |    +----------------------------------------------------------------------------------+
     |    | d4                            | Inf                                              |
     |    | d2                            | 2.0                                              |
     |    | mask                          | get_bonded_mask                                  |
     |    | d1                            | 1.0                                              |
     |    | d3                            | Inf                                              |
     |    +----------------------------------------------------------------------------------+
     |    
     └──  ○  Selection: nothing
    ```
    """
    function get_default_all_atom_clash_restraint(;α::T = ProtoSyn.Units.defaultFloat(1.0), mask::Opt{ProtoSyn.Mask} = nothing) where {T <: AbstractFloat}
        # * Note: The default :d1 and :d2 distances were parametrized based on
        # * the 2A3D PDB structure.
        
        # with the bonded mask ignore clashes between bonded atoms
        if mask === nothing
            mask = ProtoSyn.Calculators.get_bonded_mask
        end

        return EnergyFunctionComponent(
            "All_Atom_Clash_Rest",
            ProtoSyn.Calculators.Restraints.calc_flat_bottom_restraint,
            nothing,
            Dict{Symbol, Any}(:d1 => 1.0, :d2 => 2.0, :d3 => Inf, :d4 => Inf, :mask => mask),
            α,
            true)
    end
end