module Mutators

    using ProtoSyn

    abstract type AbstractMutator end

    include("dihedral.jl")
    include("crankshaft.jl")

end