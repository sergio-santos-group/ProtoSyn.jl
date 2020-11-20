module Common
    # The Common module stores function with higher levels of abstraction that
    # are common practice in ProtoSyn algorithms

    using ProtoSyn.Units
    using ProtoSyn.Calculators
    using ProtoSyn.Peptides

    export get_default_energy_function
    get_default_energy_function(::Type{T}) where {T <: AbstractFloat} = begin
        return Calculators.EnergyFunction(Dict(
            Calculators.TorchANI.torchani_model => T(1.0),
            Peptides.Calculators.Caterpillar.solvation_energy => T(0.01)
        ))
    end

    get_default_energy_function() = get_default_energy_function(Units.defaultFloat)

    # export get_energy_function
    # function get_energy_function(params::Opt{Calculators.EnergyFunctionComponents} = nothing)::Function

    #     if params === nothing
    #         params = get_default_energy_components()
    #     end

    #     return (pose::Pose) -> begin
    #         e = 0.0
    #         for (component, ɑ) in params.components
    #             if ɑ > 0.0
    #                 e_comp = ɑ * component.calc(pose)
    #                 pose.state.e[Symbol(component.name)] = e_comp
    #                 e += e_comp
    #             end
    #         end

    #         pose.state.e[:Total] = e
    #         return e
    #     end
    # end

end