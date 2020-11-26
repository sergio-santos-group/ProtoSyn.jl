module Mutators

    using ProtoSyn

    abstract type AbstractMutator end

    include("dihedral.jl")
    include("crankshaft.jl")
    include("rigid_body.jl")

end