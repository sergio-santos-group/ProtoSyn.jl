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

    if "USE_SEQ_DES" in keys(ENV) && ENV["USE_SEQ_DES"] === "false"
        @warn "Environment variable USE_SEQ_DES set to `false`. Not loading torchani."
    else
        include("seq_des.jl")
    end

    show_available_energy_function_components() = ProtoSyn.Calculators.show_available_energy_function_components(stdout, ProtoSyn.Peptides.Calculators)
    get_available_energy_function_components() = ProtoSyn.Calculators.get_available_energy_function_components(ProtoSyn.Peptides.Calculators)
end