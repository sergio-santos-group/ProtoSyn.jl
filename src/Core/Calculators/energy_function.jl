using Printf

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