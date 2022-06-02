```@meta
CurrentModule = ProtoSyn
```

# [Input and Output (IO)](@id core-io)

The following methods allow the user to load and export models and data from and to files.

By default, ProtoSyn supports the following file formats. This may change in future versions of ProtoSyn!

| Format extension | Read | Write |
| ---------------- | :--: | ----: |
| .pdb             | Yes  |   Yes |
| .yml             | Yes  |   Yes |
| .pqr             |  No  |   Yes |
| .xyz             |  No  |   Yes |

```@docs
load
download
write
append
```

The following methods are used by ProtoSyn in the [`load`](@ref) function automatically, but are displayed here for documentation purposes only.

```@docs
is_trajectory
splice_trajectory
load_trajectory
```