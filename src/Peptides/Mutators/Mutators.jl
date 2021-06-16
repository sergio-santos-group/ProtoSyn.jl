@doc """
    Mutators

Mutators are complex functions that perform 1 specific change in a system.
Besides the list of available Mutators at ProtoSyn.Mutators, the Peptides
module adds the following Mutators:
- `RotamerMutator`
- `DesignMutator`
"""
module Mutators

    using ProtoSyn
    using ProtoSyn.Mutators: AbstractMutator
    using ProtoSyn.Peptides

    include("rotamer.jl")
    include("design.jl")
end