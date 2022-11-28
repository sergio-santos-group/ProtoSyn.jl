using Printf

"""
    EnergyFunction(components::Vector{EnergyFunctionComponent}, clean_cache_every::Int16, cache::Int16, components_by_name::Dict{String, Int}, selection::AbstractSelection, update_forces::Bool)

Construct and return a new [`EnergyFunction`](@ref) instance. An
[`EnergyFunction`](@ref) is a collection of [`EnergyFunctionComponent`](@ref)
instances, where each of these components calculates an energetic contribution
to the `:total` energy and forces acting on a system. An
[`EnergyFunctionComponent`](@ref) can be retrieved by its index or by its name
(as saved in `components_by_name`). The `AbstractSelection` `selection` defines
the [`Atom`](@ref) selection this [`EnergyFunction`](@ref) is applied to. If 
inner [`EnergyFunctionComponent`](@ref) instances also have `AbstractSelection`
selections defined, the resulting selection will be the intersection between 
both. The `update_forces` sets whether to calculate and update the
[`Pose`](@ref) [State](@ref state-types) forces. The Julia cache is automatically cleaned by
garbage collection. However, in certain cases (such as using the
[TorchANI](@ref) [`EnergyFunctionComponent`](@ref)), a manual call to garbage
collection is necessary (see
[Issue 55140](https://discourse.julialang.org/t/using-gpu-via-pycall-causes-non-reusable-memory-allocation/55140/4)).
In such cases, the [`EnergyFunction`](@ref) object has an internal `cache` that
is cleaned (by calling garbage collection) at intervals of `clean_cache_every`
calls. This interval is automatically adjusted down by measuring the current
load on the GPU, calling garbage collection once the memory allocation is
greater than `ProtoSyn.Units.max_gpu_allocation`.

    EnergyFunction([::Type{T}])

Construct and return an empty [`EnergyFunction`](@ref) instance. The
`:clean_cache_every` field is set to `ProtoSyn.Units.defaultCleanCacheEvery`, 
the `AbstractSelection` `selection` field is defined as an atomic
[`TrueSelection`](@ref) and `update_forces` is set to `false.`

    EnergyFunction(components::Vector{EnergyFunctionComponent{T}}) where {T <: AbstractFloat}

Construct and return a new [`EnergyFunction`](@ref) instance filled with the
given list of [`EnergyFunctionComponent`](@ref) instances `components`. The
`:clean_cache_every` field is set to `ProtoSyn.Units.defaultCleanCacheEvery`,
the `AbstractSelection` `selection` field is defined as an atomic
[`TrueSelection`](@ref) and `update_forces` is set to `false.`

# Fields
* `components::Vector{EnergyFunctionComponent}` - The list of [`EnergyFunctionComponent`](@ref) instances in this [`EnergyFunction`](@ref);
* `clean_cache_every::Int16` - Forcefully call garbage collection every `N` calls;
* `cache::Int16` - Current number of calls performed. Resets to zero every `clean_cache_every`;
* `components_by_name::Dict{String, Int}` - The list of [`EnergyFunctionComponent`](@ref) instances in this [`EnergyFunction`](@ref), indexed by `:name`;
* `selection::AbstractSelection` - The `AbstractSelection` selecting [`Atom`](@ref) instances to apply this [`EnergyFunction`](@ref) to;
* `update_forces::Bool` - Whether to calculate and update forces with this [`EnergyFunction`](@ref).

# See also
[`EnergyFunctionComponent`](@ref)

# Examples
```
julia> energy_function = ProtoSyn.Calculators.EnergyFunction()
🗲  Energy Function (0 components):
+----------------------------------------------------------------------+
| Index | Component name                                | Weight (α)   |
+----------------------------------------------------------------------+
+----------------------------------------------------------------------+
 ● Update forces: false
 ● Selection: Set
 └── TrueSelection (Atom)

julia> push!(energy_function, Calculators.Restraints.get_default_bond_distance_restraint())
🗲  Energy Function (1 components):
+----------------------------------------------------------------------+
| Index | Component name                                | Weight (α)   |
+----------------------------------------------------------------------+
| 1     | Bond_Distance_Rest                            |      1.000   |
+----------------------------------------------------------------------+
 ● Update forces: false
 ● Selection: Set
 └── TrueSelection (Atom)

 julia> energy_function["Bond_Distance_Rest"].α = 0.5
 0.5

julia> energy_function
🗲  Energy Function (1 components):
+----------------------------------------------------------------------+
| Index | Component name                                | Weight (α)   |
+----------------------------------------------------------------------+
| 1     | Bond_Distance_Rest                            |      0.500   |
+----------------------------------------------------------------------+
 ● Update forces: false
 ● Selection: Set
 └── TrueSelection (Atom)
```
"""
mutable struct EnergyFunction
    components::Vector{EnergyFunctionComponent}
    clean_cache_every::Int16
    cache::Int16
    components_by_name::Dict{String, Int}
    selection::AbstractSelection
    update_forces::Bool
end

EnergyFunction() = begin
    return EnergyFunction(ProtoSyn.Units.defaultFloat)
end

EnergyFunction(::Type{T}) where {T <: AbstractFloat} = begin
    return EnergyFunction(Vector{EnergyFunctionComponent{T}}())
end

EnergyFunction(components::Vector{EnergyFunctionComponent{T}}) where {T <: AbstractFloat} = begin
    energy_function = EnergyFunction(Vector{EnergyFunctionComponent}(), ProtoSyn.Units.defaultCleanCacheEvery, Int16(0), Dict{String, Int}(), TrueSelection{Atom}(), false)
    for component in components
        push!(energy_function, component)
    end

    return energy_function
end

# * Overloads ------------------------------------------------------------------
Base.length(energy_function::EnergyFunction) = length(energy_function.components)

Base.push!(energy_function::EnergyFunction, component::EnergyFunctionComponent) = begin
    push!(energy_function.components, component)
    energy_function.components_by_name[component.name] = length(energy_function)
    return energy_function
end

Base.getindex(energy_function::EnergyFunction, i::Int) = energy_function.components[i]

Base.getindex(energy_function::EnergyFunction, name::String) = begin
    return energy_function.components[energy_function.components_by_name[name]]
end

Base.pop!(energy_function::EnergyFunction) = begin
    component = pop!(energy_function.components)
    delete!(energy_function.components_by_name, component.name)
    for (i, comp) in enumerate(energy_function.components)
        energy_function.components_by_name[comp.name] = i
    end
    return component
end

Base.pop!(energy_function::EnergyFunction, efc::EnergyFunctionComponent) = begin
    e = findall(x -> x == efc, energy_function.components)
    deleteat!(energy_function.components, e)
    delete!(energy_function.components_by_name, efc.name)
    for (i, comp) in enumerate(energy_function.components)
        energy_function.components_by_name[comp.name] = i
    end
    return efc
end

function Base.copy(ef::EnergyFunction)
    nef = EnergyFunction()
    nef.clean_cache_every = ef.clean_cache_every
    for component in ef.components
        push!(nef.components, copy(component))
    end
    nef.components_by_name = copy(ef.components_by_name)
    nef.selection = ef.selection
    nef.update_forces = ef.update_forces
    return nef
end


"""
    fixate_masks!(ef::EnergyFunction, pose::Pose) where {T <: AbstractFloat}

Change the current [`Mask`](@ref) type of all [`EnergyFunctionComponent`](@ref) 
instances in the given [`EnergyFunction`](@ref) `ef` from dynamic to static, by
applying them to the given [`Pose`](@ref) `pose`.

# See also
[`fixate_mask!`](@ref ProtoSyn.Calculators.fixate_mask!)

# Examples
```jldoctest
julia> ProtoSyn.Calculators.fixate_masks!(energy_function, pose)

julia> energy_function[4].settings[:mask]
ProtoSyn.Mask
 ├── Type: Atom
 ├── Size: (343, 343)
 ├── Count: 111594 / 117649
 └── Content: [0 0 … 1 1; 0 0 … 1 1; … ; 1 1 … 0 0; 1 1 … 0 0]
```
"""
function fixate_masks!(ef::EnergyFunction, pose::Pose)
    for efc in ef.components
        if (:mask in keys(efc.settings)) && isa(efc.settings[:mask], Function)
            efc.settings[:mask] = efc.settings[:mask](pose, ef.selection & efc.selection) # intersected masks
        end
    end
end

# * Call energy function -------------------------------------------------------

function (energy_function::EnergyFunction)(pose::Pose; update_forces_overwrite::Opt{Bool} = nothing)
    e              = 0.0
    performed_calc = false
    sync!(pose)

    for component in energy_function.components

        # Calculate component only if α ≠ 0.0
        if component.α !== 0.0

            # Intersect update forces
            if update_forces_overwrite !== nothing
                uf = update_forces_overwrite
            else
                uf = energy_function.update_forces & component.update_forces
            end

            # Intersect selection
            sele = energy_function.selection & component.selection

            # Reset energy & forces on first calculation. Both are reset if one
            # of them is updated.
            if !performed_calc
                pose.state.e = Dict{Symbol, eltype(pose.state)}()
                fill!(pose.state.f, eltype(pose.state)(0))
            end

            # Calculate component
            energy, forces = component.calc(
                ProtoSyn.acceleration.active, pose, sele, uf;
                component.settings...)
                
            performed_calc = true
            e_comp         = component.α * energy
            pose.state.e[Symbol(component.name)] = e_comp
            e += e_comp

            # Update forces (requires correctly indexed pose)
            if uf & !(forces === nothing)
                if sele !== nothing
                    i = 0
                    selected = sele(pose) # this is a mask
                    for atom_index in 1:pose.state.size
                        !selected[atom_index] && continue
                        i += 1
                        pose.state.f[:, atom_index] += forces[:, i] .* component.α
                    end
                else
                    for atom_index in 1:pose.state.size
                        pose.state.f[:, atom_index] += forces[:, atom_index] .* component.α
                    end
                end
            end
        end
    end

    # Perform memory allocation clean-up and maintenance
    if performed_calc && ProtoSyn.acceleration.active == ProtoSyn.CUDA_2
        energy_function.cache += 1
        if energy_function.cache % energy_function.clean_cache_every == 0
            GC.gc(false)
            energy_function.cache = 0
        else
            alloc = ProtoSyn.gpu_allocation()
            if alloc > ProtoSyn.Units.max_gpu_allocation
                GC.gc(false)
                energy_function.clean_cache_every = energy_function.cache
                energy_function.cache = 0
            end
        end
    end

    pose.state.e[:Total] = e
    return e
end

# * Show -----------------------------------------------------------------------

function Base.show(io::IO, efc::EnergyFunction, level_code::Opt{LevelCode} = nothing)
    level_code = level_code === nothing ? LevelCode() : level_code
    lead       = ProtoSyn.get_lead(level_code)
    inner_lead = ProtoSyn.get_inner_lead(level_code)

    println(io, lead*"🗲  Energy Function ($(length(efc.components)) components):")
    println(io, inner_lead*"+"*repeat("-", 70)*"+")
    @printf(io, "%s| %-5s | %-45s | %-12s |\n", inner_lead, "Index", "Component name", "Weight (α)")
    println(io, inner_lead*"+"*repeat("-", 70)*"+")
    for (index, component) in enumerate(efc.components)
        if component.α < 0.01
            @printf(io, "%s| %-5d | %-45s | %10.2e   |\n", inner_lead, index, component.name, component.α)
        else
            @printf(io, "%s| %-5d | %-45s | %10.2f   |\n", inner_lead, index, component.name, component.α)
        end
    end
    println(io, inner_lead*"+"*repeat("-", 70)*"+")

    println(io, inner_lead*" ● Update forces: $(efc.update_forces)")
    println(io, inner_lead*" ● Selection: Set")
    Base.show(io, efc.selection, vcat(level_code, 4))
end