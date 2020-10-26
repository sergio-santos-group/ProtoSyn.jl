mutable struct EnergyFunctionComponents{T <: AbstractFloat}

    components::Dict{EnergyFunctionComponent, T}
end

EnergyFunctionComponents(::Type{T}) where {T <: AbstractFloat} = begin
    EnergyFunctionComponents(Dict(
        TorchANI.torchani_model => T(1.0),
        Caterpillar.solvation_energy => T(0.01)
    ))
end

EnergyFunctionComponents() = EnergyFunctionComponents(ProtoSyn.Units.defaultFloat)

function Base.show(io::IO, efc::EnergyFunctionComponents)
    println(io, "Energy Function Components:")
    for (index, (component, ɑ)) in enumerate(efc.components)
        println(io, "$index) $(component.name) - $(ɑ)")
    end
end

function get_energy_function(params::Opt{EnergyFunctionComponents} = nothing)::Function

    if params === nothing
        params = EnergyFunctionComponents()
    end

    return (pose::Pose) -> begin
        e = 0.0
        for (component, ɑ) in params.components
            if ɑ > 0.0
                e_comp = ɑ * component.calc(pose)
                pose.state.e[Symbol(component.name)] = e_comp
                e += e_comp
            end
        end

        pose.state.e[:Total] = e
        return e
    end
end