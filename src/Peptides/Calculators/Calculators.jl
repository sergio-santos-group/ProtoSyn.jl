module Calculators

    using ProtoSyn.Peptides
    using Serialization
    using Polynomials

    phi_file = joinpath(Peptides.resource_dir, "ramachandran/phi-potential.jls")
    psi_file = joinpath(Peptides.resource_dir, "ramachandran/psi-potential.jls")
    phi_potential = open(deserialize, phi_file)
    psi_potential = open(deserialize, psi_file)

    include("Caterpillar/Caterpillar.jl")
    include("restraints.jl")
    include("hydrogen_bonds.jl")
end