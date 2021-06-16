```@meta
CurrentModule = ProtoSyn.Drivers
```

# Compound Driver

Similar to the [Compound Mutator](@ref), a [Compound Driver](@ref) joins and sequentially calls a set of `Driver` and `AbstractMutator` instances, as well as custom sampling functions (the only difference to a [`CompoundMutator`](@ref ProtoSyn.Mutators.CompoundMutator) is the addition of `Driver` instances as callable objects).

```@docs
CompoundDriver
```