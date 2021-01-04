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
end