using Printf

mutable struct EnergyFunction{T <: AbstractFloat}

    components::Dict{EnergyFunctionComponent, T}
    clean_cache_every::Int16
end

EnergyFunction(components::Dict{EnergyFunctionComponent, T}) where {T <: AbstractFloat} = begin
    return EnergyFunction{T}(components, Int16(20))
end


function (energy_function::EnergyFunction)(pose::Pose; update_forces::Bool = false)
    e = 0.0
    for (component, ɑ) in energy_function.components
        if ɑ > 0.0
            energy, forces = component.calc(pose, update_forces = update_forces)
            e_comp = ɑ * energy
            pose.state.e[Symbol(component.name)] = e_comp
            e += e_comp
            if update_forces
                for atom_index in 1:pose.state.size
                    pose.state.f[:, atom_index] += forces[atom_index, :]
                end
            end

            component.cached_n_calls += 1
            if component.cached_n_calls % energy_function.clean_cache_every == 0
                GC.gc(false)
                component.cached_n_calls = 0
            end
        end
    end

    pose.state.e[:Total] = e
    return e
end


function Base.show(io::IO, efc::EnergyFunction)
    println(io, "\n+"*repeat("-", 58)*"+")
    @printf(io, "| %-5s | %-35s | %-10s |\n", "Index", "Component name", "Weight (ɑ)")
    println(io, "+"*repeat("-", 58)*"+")
    for (index, (component, ɑ)) in enumerate(efc.components)
        @printf(io, "| %-5d | %-35s | %-10.2f |\n", index, component.name, ɑ)
    end
    println(io, "+"*repeat("-", 58)*"+")
end