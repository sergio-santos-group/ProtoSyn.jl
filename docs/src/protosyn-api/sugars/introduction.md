# Sugars

Besides the Core module, ProtoSyn comes with several other extra modules that build on top of and expand the base methods, types and submodules. These ProtoSyn modules are specific for a given topic. In the case of the [Sugars](@ref) module, all new or expanded methods and types explore the topic of carbohydrates simulation and ramified sugar structures, including glycoproteins. With this extra constraint, most methods can be especialized.

In sum, the [Sugars](@ref) module makes available several methods, types and submodules specific for ramified polysaccharide structures.

!!! ukw "Note:"
    🛠 **Under construction.** This module is not finished: several bugs and incomplete features can be found. The existing code is provided as is, and should be mainly used for further improvement and development.

In order to use this module, including the following call is often useful:

```@setup peptides
using ProtoSyn
```

```@example peptides
using ProtoSyn.Sugars
```