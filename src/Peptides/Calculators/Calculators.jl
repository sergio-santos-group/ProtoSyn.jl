module Calculators

    using ProtoSyn.Peptides
    using Serialization
    using Polynomials

    include("Caterpillar/Caterpillar.jl")
    include("restraints.jl")
    include("natural_frequency.jl")
    include("ss_propensity.jl")
    include("electrostatics.jl")
end