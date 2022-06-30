```@meta
CurrentModule = ProtoSyn
```

# REF-15

The [Rosetta](https://www.rosettacommons.org/software) ( and python wrapper [PyRosetta](https://www.pyrosetta.org/)) software has been on the forefront of technological advances in the scope of molecular simulations and energy function development in the past few decades. One such leap was the development of [REF-15](https://pubs.acs.org/doi/10.1021/acs.jctc.7b00125), an all-atom energy function for macromolecular modeling and design. ProtoSyn directly links to PyRosetta to make use of the REF-15 energy function in its simulations, as an [`EnergyFunctionComponent`](@ref ProtoSyn.Calculators.EnergyFunctionComponent).

!!! ukw "Note:"
    Given the natural differences between the way the ProtoSyn and PyRoseta packages handle a [`Pose`](@ref), currently, the only way to convert between the two data structures in through the file system. This process, as expected, carries a high computational cost. If no mutation/design process is being performed, it's possible to re-use the sabe PyRosetta pose by simply updating the atomic coordinates. Consider employing the [`fixate_rosetta_pose!`](@ref ProtoSyn.Calculators.REF15.fixate_rosetta_pose!) method.

```@docs
ProtoSyn.Calculators.REF15.calc_ref15
ProtoSyn.Calculators.REF15.get_default_ref15
ProtoSyn.Calculators.REF15.fixate_rosetta_pose!
```