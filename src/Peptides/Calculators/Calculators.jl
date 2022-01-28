module Calculators

    using ProtoSyn.Peptides
    using Serialization
    using Polynomials

    include("Caterpillar/Caterpillar.jl")
    include("restraints.jl")
    include("hydrogen_bonds.jl")
    include("natural_frequency.jl")
end