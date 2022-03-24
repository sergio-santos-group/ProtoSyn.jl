module Calculators

    using ProtoSyn.Peptides
    using Serialization
    using Polynomials

    include("Caterpillar/Caterpillar.jl")
    include("restraints.jl")
    include("natural_frequency.jl")
    include("ss_propensity.jl")
    include("electrostatics.jl")
    include("sasa.jl")

    show_available_energy_function_components() = ProtoSyn.Calculators.show_available_energy_function_components(stdout, ProtoSyn.Peptides.Calculators)
    get_available_energy_function_components() = ProtoSyn.Calculators.get_available_energy_function_components(ProtoSyn.Peptides.Calculators)
end