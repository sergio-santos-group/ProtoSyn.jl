```@meta
CurrentModule = ProtoSyn.GMX
```

# [External packages](@id external-packages)

Because ProtoSyn does not exist in an isolated bubble but, instead, in a flourishing and ever-growing sea of novel applications and tools for computational chemistry, the [External packages](@ref external-packages) submodules attempts to link structures (and simulations) from ProtoSyn to any file format or input type for common and useful external packages.

## Gromacs

[Gromacs](https://www.gromacs.org/) is "a free and open-source software suite for high-performance molecular dynamics and output analysis". ProtoSyn makes available tools to automatically export [`Pose`](@ref) structures to the correct input formats for [Gromacs](https://www.gromacs.org/) (including the definition of .itp files with the forcefield parameters for molecular dynamics), or even employ [Gromacs](https://www.gromacs.org/) tools to manipulate structures directly.

```@docs
generate_gmx_itp
add_bounding_box
add_solvent
add_ions
```