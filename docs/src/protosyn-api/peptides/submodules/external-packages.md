```@meta
CurrentModule = ProtoSyn.Peptides.GMX
```

# [External packages](@id peptides-external-packages)

> The [External packages](@ref peptides-external-packages) is a submodule of `ProtoSyn.Peptides` module. As such, the following section introduces both new types and methods that work together, in a generally independent way from the rest of the module, and require an unique exploratory section on their own.

As futher explained in the Core module [External packages](@ref external-packages) section, ProtoSyn makes available tools to define and manipulate input structures for common and useful external packages.

## Gromacs

For the Gromacs simulation suite, the Peptides module introduces methods to define .itp files for protein structures, as well as manipulate the [`Pose`](@ref) graph to match both the [`Atom`](@ref) names and [`Atom`](@ref) order expected by Gromacs. 

```@docs
generate_gmx_itp
generate_gmx_topology
sort_atoms_and_graph_gmx
GMXS
assign_gmx_atom_names!
generate_gmx_files
```