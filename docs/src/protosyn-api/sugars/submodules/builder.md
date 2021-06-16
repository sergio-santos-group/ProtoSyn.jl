```@meta
CurrentModule = ProtoSyn.Sugars
```

# [Builder](@id sugars-builder)

> The [Builder](@ref sugars-builder) is a submodule of `ProtoSyn.Sugars` module. As such, the following section introduces both new types and methods that work together, in a generally independent way from the rest of the module, and require an unique exploratory section on their own.

As an expansion of the Core module [Builder](@ref core-builder), this submodule introduces the carbohydrates [`LGrammar`](@ref) type, as well as the necessary methods to append and insert [`Fragment`](@ref) instances from a derivation.

```@docs
grammar
```

![ProtoSyn Ramified Sugar](../../../assets/ProtoSyn-ramified-sugar.png)

**Figure 1 |** Small example of a ramified sugar (such as amylopectin), built in ProtoSyn.