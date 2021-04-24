using Printf

"""
    EnergyFunction(components::Vector{EnergyFunctionComponent}, clean_cache_every::Int16, cache::Int16, components_by_name::Dict{String, Int})

Construct and return a new [`EnergyFunction`](@ref) instance. An
[`EnergyFunction`](@ref) is a collection of [`EnergyFunctionComponent`](@ref)
instances, where each of these components calculates an energetic contribution
to the `:total` energy and forces acting on a system. An
[`EnergyFunctionComponent`](@ref) can be retrieved by its index or by its name
(as saved in `components_by_name`). The Julia cache is automatically cleaned by
garbage collection. However, in certain cases (such as using the
[TorchANI](@ref) [`EnergyFunctionComponent`](@ref)), a manual call to garbage
collection is necessary (see [Issue 55140](https://discourse.julialang.org/t/using-gpu-via-pycall-causes-non-reusable-memory-allocation/55140/4)).
In such cases, the [`EnergyFunction`](@ref) object has an internal `cache` that
is cleaned (by calling garbage collection) at intervals of `clean_cache_every`
calls. This interval is automatically adjusted down by measuring the current
load on the GPU, calling garbage collection once the memory allocation is
greater than `ProtoSyn.Units.max_gpu_allocation`.

    EnergyFunction([::Type{T}])

Construct and return an empty [`EnergyFunction`](@ref) instance. The
`:clean_cache_every` field is set to `ProtoSyn.Units.defaultCleanCacheEvery`.

    EnergyFunction(components::Vector{EnergyFunctionComponent{T}}) where {T <: AbstractFloat}

Construct and return a new [`EnergyFunction`](@ref) instance filled with the
given list of [`EnergyFunctionComponent`](@ref) instances `components`. The
`:clean_cache_every` field is set to `ProtoSyn.Units.defaultCleanCacheEvery`.  

# Fields
* `components::Vector{EnergyFunctionComponent}` - The list of [`EnergyFunctionComponent`](@ref) instances in this [`EnergyFunction`](@ref);
* `clean_cache_every::Int16` - Forcefully call garbage collection every `N` calls;
* `cache::Int16` - Current number of calls performed. Resets to zero every `clean_cache_every`;
* `components_by_name::Dict{String, Int}` - The list of [`EnergyFunctionComponent`](@ref) instances in this [`EnergyFunction`](@ref), indexed by `:name`.

# See also
[`EnergyFunctionComponent`](@ref)

# Examples
```jldoctest
julia> energy_function = Calculators.EnergyFunction()

+----------------------------------------------------------+
| Index | Component name                      | Weight (α) |
+----------------------------------------------------------+
+----------------------------------------------------------+

julia> push!(energy_function, Calculators.Restraints.get_default_bond_distance_restraint())

+----------------------------------------------------------+
| Index | Component name                      | Weight (α) |
+----------------------------------------------------------+
| 1     | Bond_Distance_Restraint             | 1.000      |
+----------------------------------------------------------+

julia> energy_function["Bond_Distance_Restraint"].α = 0.5
0.5

julia> energy_function

+----------------------------------------------------------+
| Index | Component name                      | Weight (α) |
+----------------------------------------------------------+
| 1     | Bond_Distance_Restraint             | 0.500      |
+----------------------------------------------------------+
```
"""
mutable struct EnergyFunction

    components::Vector{EnergyFunctionComponent}
    clean_cache_every::Int16
    cache::Int16
    components_by_name::Dict{String, Int}
end

EnergyFunction() = begin
    return EnergyFunction(ProtoSyn.Units.defaultFloat)
end

EnergyFunction(::Type{T}) where {T <: AbstractFloat} = begin
    return EnergyFunction(Vector{EnergyFunctionComponent{T}}())
end

EnergyFunction(components::Vector{EnergyFunctionComponent{T}}) where {T <: AbstractFloat} = begin
    energy_function = EnergyFunction(Vector{EnergyFunctionComponent}(), ProtoSyn.Units.defaultCleanCacheEvery, Int16(0), Dict{String, Int}())
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
    return component
end

function Base.copy(ef::EnergyFunction)
    nef = EnergyFunction()
    nef.clean_cache_every = ef.clean_cache_every
    for (component, α) in ef.components
        nef.components[component] = α
    end
    return nef
end

# * Call energy function -------------------------------------------------------

function (energy_function::EnergyFunction)(pose::Pose, update_forces::Bool = false)
    e = 0.0
    performed_calc = false
    for component in energy_function.components
        uf = update_forces & component.update_forces
        if component.α > 0.0
            if !performed_calc 
                pose.state.e = Dict{Symbol, eltype(pose.state)}()
                uf && fill!(pose.state.f, eltype(pose.state)(0))
            end
            energy, forces = component.calc(pose, uf; component.settings...)
            performed_calc = true
            e_comp         = component.α * energy
            pose.state.e[Symbol(component.name)] = e_comp
            e += e_comp

            if uf & !(forces === nothing)
                for atom_index in 1:pose.state.size
                    pose.state.f[:, atom_index] += forces[:, atom_index] .* component.α
                end
            end
        end
    end

    # Perform memory allocation clean-up and maintenance
    if performed_calc
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

function Base.show(io::IO, efc::EnergyFunction)
    println(io, "\n+"*repeat("-", 58)*"+")
    @printf(io, "| %-5s | %-35s | %-10s |\n", "Index", "Component name", "Weight (α)")
    println(io, "+"*repeat("-", 58)*"+")
    for (index, component) in enumerate(efc.components)
        @printf(io, "| %-5d | %-35s | %-10.3f |\n", index, component.name, component.α)
    end
    println(io, "+"*repeat("-", 58)*"+")
end