module RG

    using ProtoSyn

    """
    calc_radius_gyration(pose::Pose, [selection::Opt{AbstractSelection} = nothing])

    Calculate each dimension (X, Y and Z) radius of gyration (the root mean
    square deviation to all of the structure's [`Atom`](@ref) instances from the
    [`Pose`](@ref) `pose` [`center_of_mass`](@ref ProtoSyn.center_of_mass)).
    Each [`Atom`](@ref) type mass is retrieved from `ProtoSyn.Units.mass`.

    # See also
    [`calc_radius_gyration_energy`](@ref) [`get_default_rg`](@ref)

    # Examples
    ```jldoctest
    julia> ProtoSyn.Calculators.SASA.calc_radius_gyration(pose)
    3×1 Matrix{Float64}:
    4.976766149665673
    5.501503499039599
    10.74091335737219
    ```
    """
    function calc_radius_gyration(pose::Pose, selection::Opt{AbstractSelection} = nothing)
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
    end # function


    """
        calc_radius_gyration_energy(pose::Pose, selection::Opt{AbstractSelection}, update_forces::Bool)

    Calculate the sum of all dimensions (X, Y and Z) radius of gyration. This
    Calculator does not calculate forces. As such, `update_forces` has no effect
    and exists only in order to standardize calls between Calculators.

    # See also
    [`calc_radius_gyration`](@ref) [`get_default_rg`](@ref)

    # Examples
    ```jldoctest
    julia> ProtoSyn.Calculators.SASA.calc_radius_gyration_energy(pose, nothing, false)
    (21.21918300607746, nothing)
    ```
    """
    function calc_radius_gyration_energy(::Type{A}, pose::Pose, selection::Opt{AbstractSelection}, update_forces::Bool) where {A <: ProtoSyn.AbstractAccelerationType}
        rg = calc_radius_gyration(pose, selection)
        return sum(rg), nothing
    end # function


    """
        get_default_rg(;[α::T = 1.0]) where {T <: AbstractFloat}

    Return the default Radius of Gyration energy
    [`EnergyFunctionComponent`](@ref ProtoSyn.Calculators.EnergyFunctionComponent). `α` sets
    the component weight (on an
    [`EnergyFunction`](@ref ProtoSyn.Calculators.EnergyFunction) instance, `1.0`
    by default). This function employs [`calc_radius_gyration_energy`](@ref) as
    the `:calc` function.

    # See also
    [`calc_radius_gyration_energy`](@ref)

    # Examples
    ```jldoctest
    julia> ProtoSyn.Calculators.SASA.get_default_sasa_energy()
    🞧  Energy Function Component:
    +---------------------------------------------------+
    | Name           | Radius_Gyration                  |
    | Alpha (α)      | 1.0                              |
    | Update forces  | false                            |
    | Calculator     | calc_radius_gyration_energy      |
    +---------------------------------------------------+
    └──  ○  Selection: nothing
    ```
    """
    function get_default_rg(;α::T = 1.0) where {T <: AbstractFloat}
        return ProtoSyn.Calculators.EnergyFunctionComponent(
            "Radius_Gyration",
            calc_radius_gyration_energy,
            nothing,
            Dict{Symbol, Any}(),
            α,
            false)
    end # function

end # module