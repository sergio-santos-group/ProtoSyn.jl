```@meta
CurrentModule = ProtoSyn.Peptides
```

# [Input and Output (IO)](@id peptides-io)

The following methods allow the user to load and export models and data from and to files, and build on top of the Core methods (See [Input and Output (IO)](@ref core-io)). The main different in comparison with the Core module IO methods is the addition of parenthood inference and graph sorting according to IUPAC conventions.

```@docs
load
download
```