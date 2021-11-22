@doc """
    Mutators

Mutators are complex functions that perform 1 specific change in a system.
Currently, ProtoSyn-base makes available the following mutators:
- `DihedralMutator`
- `CrankshaftMutator`
- `RotationRigidBodyMutator`
- `TranslationRigidBodyMutator`
- `BackrubMutator`
- `CompoundMutator`
"""
module Mutators

    using ProtoSyn

    abstract type AbstractMutator end

    include("dihedral.jl")
    include("crankshaft.jl")
    include("rigid_body.jl")
    include("backrub.jl")
    include("compound.jl")
end