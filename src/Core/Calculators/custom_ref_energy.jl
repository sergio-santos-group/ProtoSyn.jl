"""
    Calculators.calc_custom_ref_energy([::A], pose::Pose, selection::Opt{AbstractSelection}, update_forces::Bool = false; [map::Opt{Dict{AbstractSelection, T}} = nothing]) where {A, T <: AbstractFloat}
    
Calculate and return the [`Pose`](@ref) `pose` energy according to the given
`map`. The `map` is a 1 to 1 correspondence between `AbstractSelection`
instances and energy values. If provided, an `AbstractSelection` `selection`
limits the selected range of [`Atom`](@ref) instances considered for this
[`EnergyFunctionComponent`](@ref) calculation. Note that any provided
`selection` and `map` `AbstractSelection` type are promoted to the biggest
value: `selection` is promoted to the `AbstractSelection` type of `map` if the
`map` type is larger or vice-versa. `update_forces` `Bool` has no effect in this
[`EnergyFunctionComponent`](@ref) calculation and is only included to provide
a standard function signature between all [`EnergyFunctionComponent`](@ref)
instances. If no `map` is provided, returns 0.0.

# Examples
```
julia> ProtoSyn.Calculators.calc_custom_ref_energy(pose, nothing, false)
(0.0, nothing)

julia> ProtoSyn.Calculators.calc_custom_ref_energy(pose, nothing, false, map = Dict{AbstractSelection, Float64}(rn"ALA" => 3.0, rn"GLU" => 2.4))
(71.4, nothing)
```
"""
function calc_custom_ref_energy(A::Type{<: ProtoSyn.AbstractAccelerationType}, pose::Pose, selection::AbstractSelection, update_forces::Bool = false; map::Opt{Dict{AbstractSelection, T}} = nothing) where {T <: AbstractFloat}

    map === nothing && return 0.0, nothing
    selection_type = ProtoSyn.selection_type(selection)

    e = 0.0
    for (map_sele, weight) in map
        map_sele_type = ProtoSyn.selection_type(map_sele)
        if map_sele_type > selection_type
            sele = ProtoSyn.promote(selection, map_sele_type) & map_sele
        elseif selection_type > map_sele_type
            sele = ProtoSyn.promote(map_sele, selection_type) & selection
        end

        e += count(sele(pose).content) * weight
    end

    return e, nothing
end

function calc_custom_ref_energy(A::Type{<: ProtoSyn.AbstractAccelerationType}, pose::Pose, selection::Nothing, update_forces::Bool = false; map::Opt{Dict{AbstractSelection, T}} = nothing) where {T <: AbstractFloat}
    calc_custom_ref_energy(A, pose, TrueSelection{Atom}(), update_forces, map = map)
end

calc_custom_ref_energy(pose::Pose, selection::Opt{AbstractSelection}, update_forces::Bool = false; map::Opt{Dict{AbstractSelection, T}} = nothing) where {T <: AbstractFloat} = begin
    calc_custom_ref_energy(ProtoSyn.acceleration.active, pose, selection, update_forces, map = map)
end


"""
    get_default_custom_ref_energy(;[α::T = 1.0]) where {T <: AbstractFloat}

Return the default custom reference energy [`EnergyFunctionComponent`](@ref).
`α` sets the component weight (on an
[`EnergyFunction`](@ref ProtoSyn.Calculators.EnergyFunction) instance). This
component employs the [`calc_custom_ref_energy`](@ref) method, therefore
defining a [`Pose`](@ref) energy based on a user-defined map between
`AbstractSelection` instances and energy values.

# Settings
* `map::Dict{AbstractSelection, T}` - Defines which map the custom reference energy should use;

# Examples
```jldoctest
julia> ProtoSyn.Calculators.get_default_custom_ref_energy()
🞧  Energy Function Component:
+---------------------------------------------------+
| Name           | Custom_Ref_Energy                |
| Alpha (α)      | 1.0                              |
| Update forces  | true                             |
| Calculator     | calc_custom_ref_energy           |
+---------------------------------------------------+
 |    +----------------------------------------------------------------------------------+
 ├──  ● Settings                      | Value                                            |
 |    +----------------------------------------------------------------------------------+
 |    | map                           | Dict{AbstractSelection, Float64}(0 components)   |
 |    +----------------------------------------------------------------------------------+
 |    
 └──  ○  Selection: nothing
```
"""
function get_default_custom_ref_energy(;α::T = 1.0) where {T <: AbstractFloat}
    return EnergyFunctionComponent(
        "Custom_Ref_Energy",
        calc_custom_ref_energy,
        nothing,
        Dict{Symbol, Any}(:map => Dict{AbstractSelection, T}()),
        α,
        true)
end