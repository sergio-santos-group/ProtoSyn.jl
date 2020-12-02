module Mutators

    using ProtoSyn
    using ProtoSyn.Mutators: AbstractMutator
    using ProtoSyn.Peptides

    include("rotamer.jl")
    include("design.jl")

end