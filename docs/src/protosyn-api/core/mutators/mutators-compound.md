```@meta
CurrentModule = ProtoSyn.Mutators
```

# Compound Mutator

The [Compound Mutator](@ref) is an aggregating type of `AbstractMutator`. The sole function of this mutator is to apply multiple inner `AbstractMutator` instances in sequence. Optionally, if a selection is provided, all inner `AbstractMutator` instance's selection is merged with this common selection.

```@docs
CompoundMutator
```

![ProtoSyn Compound Mutator](../../../assets/ProtoSyn-compound-mutator.png)

**Figure 1 |** A diagram of an example [`CompoundMutator`](@ref) instance. In this example, the [`CompoundMutator`](@ref)`.mutators` list is comprised of a [`RotationRigidBodyMutator`](@ref) and a [`TranslationRigidBodyMutator`](@ref), both without any internal selections. The [`CompoundMutator`](@ref) itself, however, has a [`RandomSelectionFromList`](@ref ProtoSyn.RandomSelectionFromList) `AbstractSelection` which will randomly select one of the two provided [`RangeSelection`](@ref) instances and apply both the rotation and translation rigid body mutators.