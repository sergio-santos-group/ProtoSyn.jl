module Mutators

    using ProtoSyn

    abstract type AbstractMutator end

    include("dihedral.jl")

end