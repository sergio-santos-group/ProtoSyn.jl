![logo](./assets/logo-white.png)

# Welcome to the ProtoSyn.jl documentation!
_[ProtoSyn.jl](https://github.com/sergio-santos-group/ProtoSyn.jl) is a platform
for molecular manipulation and simulation, with an emphasis on peptide design._


The main goal of ProtoSyn is to be a basis on top of which new tools and
protocols can be experimented and prototyped. Taking advantage of Julia’s
environment, ProtoSyn has been built with emergent technologies in mind, such as
distributed computing, GPU and SIMD acceleration and machine learning models
usage.

### Features
* Edit a peptide structure by removing, adding and mutating any number of residues
* Create peptides from scratch by providing the desired sequence
* Copy parts of other molecules and graft them together to create something new
* Run Monte Carlo simulations to optimize a structure or a sequence
* Calculate energies and forces using [TorchANI](https://github.com/aiqm/torchani) code, integrated in Julia.
* Perform Steepest Descent optimizations
* Explore rotamer libraries to optimize sidechain packaging
* Select residues based on name, index, distance and other parameters, with a rich combinatory selection syntax
* Perform rigid body docking of ligands

Version 1.0 is scheduled to be released in the end of 2021.

## Getting started
* [Installation](@ref)
* [First steps](@ref)

## Examples
* Example 1
* Example 2

# ProtoSyn.jl API
