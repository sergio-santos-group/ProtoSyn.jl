# Examples

A major focus of ProtoSyn was the compilation of an extensive list of tutorials and examples. The main list of example notebook files can be found at the [Examples folder](https://github.com/sergio-santos-group/ProtoSyn.jl/tree/master/examples), in the GitHub repository. In the next section, a quick description and link to each individual example can be found.

+ [1 - Getting started](@ref examples-getting-started)
+ [2 - Selections](@ref)
+ [3 - Molecular manipulation](@ref)
+ [4 - Energy calculation](@ref)
+ [5 - Monte Carlo](@ref)
+ [6 - Steepest Descent](@ref)
+ [7 - ILS](@ref)
+ [8 - Design](@ref)
+ [9 - Sidechain packaging](@ref)
+ [10 - Rigid body docking](@ref)
+ [11 - Distributed computing](@ref)
+ [12 - Ramified polymers & sugars](@ref)
+ [13 - Advanced mutators](@ref)
+ [14 - Non-canonical aminoacids (NCAAs)](@ref)
+ [15 - Post-translational modifications](@ref)
+ [16 - Materials & carbons](@ref)

## [1 - Getting started](@id examples-getting-started)

[You can find the example file here](https://github.com/sergio-santos-group/ProtoSyn.jl/blob/master/examples/1-getting-started.ipynb).
In this first example script, we will explore how ProtoSyn.jl is organized and what are the available data structures. The example script is divided in 4 parts: loading a PDB file, exploring the graph structure, exploring the state structure and exporting a [`Pose`](@ref) into a PDB file. For more details, check the Core [Types](@ref core-types) section of this manual.

## 2 - Selections

[You can find the example file here](https://github.com/sergio-santos-group/ProtoSyn.jl/blob/master/examples/2-selections.ipynb).
ProtoSyn comes equipped with a powerful selecting syntax, useful for highlighting and specifying targets for the several manipulation tools of ProtoSyn. In this example, we will take a closer look at the different type of selections and how to apply them. For more information, check the [Selections](@ref core-selections) section of this manual.

## 3 - Molecular manipulation

[You can find the example file here](https://github.com/sergio-santos-group/ProtoSyn.jl/blob/master/examples/3-molecular-manipulation.ipynb).
One of the main applications of ProtoSyn is the easy manipulation of molecular structures: add and remove residues and whole loops, mutate aminoacids, bond and unbond atoms or completly remove them from the structure, and even create completly new sequences of aminoacids from scratch. These are just some examples of the manipulations possible with the ProtoSyn framework. In the next examples we will explore a few of these tasks. For more details, check the Core [Pose](@ref pose-methods) methods and the Peptides [Pose](@ref peptides-pose-methods) methods sections.

## 4 - Energy calculation

[You can find the example file here](https://github.com/sergio-santos-group/ProtoSyn.jl/blob/master/examples/4-energy-calculation.ipynb).
A crucial feature of the ProtoSyn framework is the ability to evaluate a given Pose based on a multitude of factors and potentials. In this example, a closer look on how to build an energy function is provided. Some of the available energy function components include Coulomb potentials, SASA/GB models, distance-based potentials, hydrogen network evaluators, etc. For more details and an exhaustive list of available potentials, check the [Calculators section](@ref).

## 5 - Monte Carlo

[You can find the example file here](https://github.com/sergio-santos-group/ProtoSyn.jl/blob/master/examples/5-monte-carlo.ipynb).
Besides the possibility to create custom simulation algorithms, ProtoSyn already provides several of these algorithms by default. A common simulation type is the Monte Carlo sampling. In this type of simulation, the conformational space is explored (or sampled) and the resulting candidate structure is evaluated and compared to the previous state. If the new energy is lower than the previous value (or accepted by an heuristic criterium, such as the Metropolis criterium), the new structure is saved for the continuation of the simulation. Otherwise, recover the previous structure and attempt a diferent change. This loop continues for N steps until the completion of the simulation. In this example, we will build a Monte Carlo simulation driver, while taking a first look at some important ProtoSyn components: the Mutators, Drivers and Callback instances. For more information, check the following sections of the manual: the [Mutators section](@ref), the [Monte Carlo Driver](@ref) section and [Callbacks](@ref) section.

## 6 - Steepest Descent

[You can find the example file here](https://github.com/sergio-santos-group/ProtoSyn.jl/blob/master/examples/6-steepest-descent.ipynb).
ProtoSyn makes available the Steepest Descent Driver, whose simulation algorithm calculates the forces being felt on each atom of a molecular structure (via an Energy Function instance) and updates the atoms position in accordance, as to relax the structure. For more information check the [Steepest Descent Driver](@ref) section.

## 7 - ILS

[You can find the example file here](https://github.com/sergio-santos-group/ProtoSyn.jl/blob/master/examples/7-ils.ipynb).
ILS stands for "Iterated Local Search" and it's a conformational search algorithm. In short, this algorithm is divided in two main components, the inner and outer loop. The inner loop can be, for example, a Monte Carlo simulation, and should explore and sample the local minimum, while the outer loop performs a "jump": a large conformational change, taking the system to a whole different local minimum for the next inner loop iteration. In this way, the conformational space can be efficiently sampled, with a higher certainty of exploring the global minimum during the simulation. In ProtoSyn, this algorithm is employed as a Driver, and will be briefly explored in this example. For more details, check the [ILS Driver](@ref) of this manual.

## 8 - Design

[You can find the example file here](https://github.com/sergio-santos-group/ProtoSyn.jl/blob/master/examples/8-design.ipynb).
One of the main goals of ProtoSyn is to provide an easy interface for the design of small peptides. This is achieved by mutating residues (more often than not in a specific region, i.e.: an active site) in order to stabilize certain interactions. In this example we will explore the DesignMutator (which provides random mutations in a selection). More more details, check the [Design Mutator](@ref) section of the manual.

## 9 - Sidechain packaging

[You can find the example file here](https://github.com/sergio-santos-group/ProtoSyn.jl/blob/master/examples/9-sidechain-packaging.ipynb).
A sub-problem of protein design is the correct packaging of aminoacid sidechains, that is, to find the correct rotamer (i.e.: set of sidechain chi dihedral angles) that minimizes clashes and augments stabilizing interactions with other sidechains in the 3D space. The conformational space to explore is, therefore, enormous. Several rotamer libraries have been proposed in the past to minimize this space by imposing certain restrictions. By default, ProtoSyn employs the Dunbrack Rotamer Library 2011, which reduces the rotamers to the most observed combination in natural databases, as well as imposing backbone dependency (i.e.: certain rotamers are only present for a given combination of backbone phi and psi dihedral angles). This greatly reduces the conformational space to search, while improving the likelihood of acceptance of a new rotamer. With this mind, in ProtoSyn, we can load a rotamer library and sample new rotamers, in a Monte Carlo simulation, in order to improve the sidechain packaging of a peptide. For more details, check the [Rotamer Mutator](@ref) section of this manual.

## 10 - Rigid body docking

[You can find the example file here](https://github.com/sergio-santos-group/ProtoSyn.jl/blob/master/examples/10-rigid-body-docking.ipynb).
ProtoSyn makes available two Rigid Body Mutators: translational and rotational. These allow the translation and rotation of blocks of residues (such as a ligand or a second molecule). With these, users are able to create custom Rigid Body Docking algorithms (for example, by moving a ligand to a grid of starting positions). In this quick example, we will take a look at a simplistic approach to rigid body docking, using a Monte Carlo simulation. In each step, short rotation and translation movements are attempted and, after evaluation by an energy function, accepted or rejected based on the Metropolis Criterium. For more information, consider reading the [Rigid Body Mutators](@ref) section of this manual.

## 11 - Distributed computing

[You can find the example file here](https://github.com/sergio-santos-group/ProtoSyn.jl/blob/master/examples/11-distributed-computing.ipynb).
By being developed in Julia, ProtoSyn enjoys some of the features naturally provided by the language, such as easy SIMD and GPU acceleration, a rich package environment and, among others, access to high level parallel and distributed computing routines. In this example we will take a look on how to launch and gather several decoys of the same Monte Carlo simulation. By having multiple decoys in parallel, given the random nature of the algorithm, we can have a greater confidence in the complete sampling of the conformational space, and that the obtained result is real. For more information, reffer to the Julia's [Distributed Computing](https://docs.julialang.org/en/v1/stdlib/Distributed/) manual section.

## 12 - Ramified polymers & sugars

[You can find the example file here](https://github.com/sergio-santos-group/ProtoSyn.jl/blob/master/examples/12-materials-sugars.ipynb).
Besides the Core and Peptides modules, ProtoSyn includes extra modules, such as the Sugars modules. In this example, we will briefly explore the goals of this modules: the introduction of ramified polymers. Note that this module is not the focus of ProtoSyn, in its current version, and may be improved in future iterations. For more details, check the [Sugars](@ref) sections of this manual.

## 13 - Advanced mutators

[You can find the example file here](https://github.com/sergio-santos-group/ProtoSyn.jl/blob/master/examples/13-advanced-mutators.ipynb).
Using the base ProtoSyn methods and mutators, users can build arbitrarily advanced and complex mutators for specific applications. In this example, we will explore how to build a Block Rotation Mutator (also known as _blockrot_). This advanced mutator is comprised of two rigid body mutators (translational and rotational), acting on a random block of sequential aminoacids, selected from a list of possible blocks. After rotation, since bonds are broken, a quick steepest descent drivers run is performed to re-close the severed connections. For this reason, the applied rotational and translational movement needs to be small. The final result is a concerted rotation of a whole block of aminoacids, in a peptide, whose occurrence using simple dihedral rotations was virtually impossible.

## 14 - Non-canonical aminoacids (NCAAs)

[You can find the example file here](https://github.com/sergio-santos-group/ProtoSyn.jl/blob/master/examples/14-non-canonical-aminoacids.ipynb).
ProtoSyn.jl allows users to define non-canonical aminoacids (NCAAs) and include them in simulations, such as design efforts and mutations. Besides using the NCAAs included with ProtoSyn, by default, in this brief tutorial the process to add a new NCAA is explored, as well as its implementation in a simple mutation application. The example NCAA is retrieved directly from a protein pose, by extracting the residue of interest from a downloaded file from Protein Data Bank (PDB).

## 15 - Post-translational modifications

[You can find the example file here](https://github.com/sergio-santos-group/ProtoSyn.jl/blob/master/examples/15-post-translational-modifications.ipynb). Post translational modifications (PTMs) are structural modifications to a structure, usually performed after synthesis, and can include the removal or addition of functional groups. Some examples of PTMs include methylation and phosphorylation processes, where methyl- and phosphoryl- groups are grafted to a given residue. In this tutorial, both these processes will be exemplified, making use of ProtoSyn's structural manipulation methods.


## 16 - Materials & carbons

[You can find the example file here](https://github.com/sergio-santos-group/ProtoSyn.jl/blob/master/examples/16-materials-carbons.ipynb). 
The objective of this example is to introduce methods and types for the creation and manipulation of both cristaline & amorphous inorganic structures. In the current version, this is restricted to a few basic unit cell repetitions and carbon sheets. These are, however, useful when setting up benchmarks for the measurement of speed/accuracy of new or improved energy functions, for example. In this example, we will create a repetition of one such basic unit cell, the primitive lattice. More interestingly, perhaps, are the methods for carbon sheet generation and manipulation. In this example, the topics of generation of simple carbon sheets, multi-layer carbon microcrystallites, pore generation and functionalization are explored and showcased.

!!! ukw "Note:"
    Would you like to see a topic covered in a more in-depth tutorial? Don't hesitate to contact or drop an issue in the [project's GitHub page](https://github.com/sergio-santos-group/ProtoSyn.jl/issues).