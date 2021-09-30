"""
    EnergyFunctionComponent(name::String, calc::Function, settings::Dict{Symbol, Any}, α::T, update_forces::Bool)

Return a new [`EnergyFunctionComponent`](@ref) instance with the given `name`.
An [`EnergyFunctionComponent`](@ref) object is responsible to calculate an
energetic contribution to the system, based on a given interaction, model or
restraint (for example, a *Bond Distance Restraint* is responsible to calculate
the energetic contribution by bonds with lengths superior to a given threshold).
The `:total` energy of a [`Pose`](@ref) is, therefore, the sum of all
[`EnergyFunctionComponent`](@ref) applied in a given calculation. Each 
[`EnergyFunctionComponent`](@ref) has a `calc` `Function` that receives a
[`Pose`](@ref) and outputs an energy value. Optionally, this `calc` `Function`
can also return a set of forces felt on all considered atoms, based on the
potential or model used, which is calculated if the `update_forces` flag is set
to `true`. For a more detailed explanation see the
[Creating a custom EnergyFunctionComponent](@ref) section of the documentation).
Additionally, certain [`EnergyFunctionComponent`](@ref) instances can have
a `settings` dictionary, containing kwargs that are passed to the `calc`
`Function`. These usually parametrize and fine tune the calculations perform
(such as setting the flat bottom restraint distances, as an example). When
calling an [`EnergyFunctionComponent`](@ref) calculation from an
[`EnergyFunction`](@ref), the relative weight of this
[`EnergyFunctionComponent`](@ref) in relation to others is given by the `α`
field.

# Fields
* `name::String` - The name of this [`EnergyFunctionComponent`](@ref). Used to index by `name` when in an [`EnergyFunction`](@ref);
* `calc::Function` - The calculation `Function` used to calculate this contribution;
* `settings::Dict{Symbol, Any}` - (Optional) A dictionary of kwargs provided to the `calc` `Function`, parameterizing its usage;
* `α::T` - The relative weight of this [`EnergyFunctionComponent`](@ref) when in an [`EnergyFunction`](@ref);
* `update_forces::Bool` - Toggle forces calculation by this [`EnergyFunctionComponent`](@ref) when in an [`EnergyFunction`](@ref). 

# See also
[`EnergyFunction`](@ref)

# Examples
```
julia> ProtoSyn.Calculators.Restraints.get_default_bond_distance_restraint()
         Name : Bond_Distance_Restraint
   Weight (α) : 1.0
Update forces : true
      Setings :
          :x0 => 2.0
```
"""
mutable struct EnergyFunctionComponent{T <: AbstractFloat}

    name::String
    calc::Function
    settings::Dict{Symbol, Any}
    α::T
    update_forces::Bool
end

function Base.copy(efc::EnergyFunctionComponent{T}) where {T <: AbstractFloat}
    return EnergyFunctionComponent{T}(
        efc.name,
        efc.calc,
        copy(efc.settings),
        copy(efc.α),
        copy(efc.update_forces))
end

"""
    # TODO
"""
function fixate_mask!(efc::EnergyFunctionComponent{T}, pose::Pose) where {T <: AbstractFloat}
    if (:mask in keys(efc.settings)) & isa(efc.settings[:mask], Function)
        efc.settings[:mask] = efc.settings[:mask](pose)
    end
end

# * Show -----------------------------------------------------------------------

function Base.show(io::IO, efc::EnergyFunctionComponent{T}) where {T <: AbstractFloat}
    @printf(io, "%14s : %-s\n", "Name", efc.name)
    @printf(io, "%14s : %-s\n", "Weight (α)", efc.α)
    @printf(io, "%14s : %-s\n", "Update forces", string(efc.update_forces))
    if length(efc.settings) == 0
        @printf(io, "%14s : -", "Setings")
    else
        @printf(io, "%14s :\n", "Setings")
        for (key, value) in efc.settings
            @printf(io, "%15s => ", ":"*string(key))
            if typeof(value) == Matrix{T}
                print(io, "Matrix{$T}($(size(value))\n")
            else
                print(io, "$value\n")
            end
        end
    end
end