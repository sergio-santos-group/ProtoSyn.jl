```@meta
CurrentModule = ProtoSyn.Mutators
```

# Mutators section

Each module in ProtoSyn may include a [Mutators section](@ref), containing types and methods allowing ProtoSyn to change a [`Pose`](@ref) instance's from one [`State`](@ref) to another. In other words, a Mutator (of type `AbstractMutator`) is responsible for implementing a given conformational change in a [`Pose`](@ref) (for example, implementing [Dihedral rotations](@ref), among others). In the context of a simulation, this is usually refered as "sampling".


![ProtoSyn Mutators](../../../assets/ProtoSyn-mutators.png)

**Figure 1 |** A diagram of an example [`CompoundDriver`](@ref), bringing together multiple `AbstractMutator` instances. In this example, both the [`DihedralMutator`](@ref) and the [`CrankshaftMutator`](@ref) are represented. Each `AbstractMutator` can be parametrized by a collection of settings, explained in each Mutator's page in this manual.

# Creating custom mutators

Mutators are most useful when employed in a [`Driver`](@ref) instance, who expects the following signature. Note that the Mutator instances are `struct` instances (of type `<: AbstractMutator`), but are also _functors_, meaning that are callable by the encompassing [`Driver`](@ref). For this call, only a [`Pose`](@ref) `pose` instance is given as input, and all conformational changes should be in-place (no output is expected).

```
struct NewMutator <: AbstractMutator
    # content
end

function (mutator::NewMutator)(pose::Pose)
    # content
end
```

# Available mutators

Besides all the machinery to create custom mutators, ProtoSyn's `Core` module makes available some simple `AbstractMutator` instances. Since this is not a specialized module, available mutators should be indiferent to the type of molecular structure being sampled. As of ProtoSyn 1.0, the following `AbstractMutator` instances are made available from the `Core` module, and can be further explored in the respectivelly dedicated pages of this manual:

* [Dihedral Mutator](@ref)