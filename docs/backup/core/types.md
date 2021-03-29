```@meta
CurrentModule = ProtoSyn
```

# Types

This section list all core types defined in ProtoSyn.

```@contents
```

!!! note

    This is a Note


## Containers

All container types are subtypes of AbstractContainer{T} and AbstractDigraph

```@setup abc
using ProtoSyn
```

```@repl abc
Atom <: ProtoSyn.AbstractContainer
Atom <: ProtoSyn.AbstractDigraph
```

```@docs
Atom
```

## Abstract types

```@docs
AbstractContainer
```
