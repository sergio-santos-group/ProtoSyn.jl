# First steps

## Introduction

The rational placement of amino acids in a sequence directly correlates to the 3D structure of the peptide once folded. This folded structure, in turn, dictates the interactions with the environment and therefore the function of the peptide. Being able to design new peptides for specific functions would unlock the potential of unknown conformations not yet explored by nature, with applications in medical fields, agriculture, biological remediation, enzymatic synthesis, among others. 

This was traditionally explored by random blind mutagenesis which is an expensive and time intensive experimental practice. With the evolution of computational power over the last couple of decades, computational design of small proteins has become the focus of scientific breakthroughs. Using computer simulations saves precious time and monetary costs of experiments, focusing efforts on simulated prototypes that have shown promising results. As such, multiple software solutions have been proposed over the years. An example would be [Rosetta](https://www.rosettacommons.org/software) (and its Python wrapper, [PyRosetta](http://www.pyrosetta.org/)), which has been indisputably invaluable as a platform for molecular manipulation and peptide design.

However, as it happens with so many scientific-purposed software packages, Rosetta has fallen into the two-language problem. The core of the simulation code is written in C with a more user-friendly wrap in Python that exposes some of the functionality. This has gravelly impaired the community ability to upgrade and modify this package, as well as imposing a steep learning curve to non-specialized would-be-users.

ProtoSyn, empowered by the Julia language ecosystem, aims to put forward a simple and easy to use platform for molecular manipulation and peptide design. A Julia-based solution to this challenge would naturally benefit from the native features of the language, such as easy parallelization and distributed computing, GPU acceleration and machine learning tools, among others.
  

## How is ProtoSyn organized?

![organization](../assets/ProtoSyn-organization.png)

ProtoSyn's main `struct` is the [Pose], which holds all the required information regarding a molecular system. This information is divided into a **Directional graph** and a **State**.

A directional graph is simply an hierarchical structuration of the molecular system:
At the top level, a [Topology] can hold multiple [Segments] \(which are contiguous chains of a molecule), which in turn can hold multiple residues (such as amino acids, in the case of proteins and peptides), which, finally, hold one or multiple [Atom] instances.

Each of these levels are called [AbstractContainers] and can have defining parameters. In the case of atoms, for example, an [Atom] is described by an `:id`, `:symbol`, `:name`, etc.

This graph is called **directional** because each of these components has a _parent_ and can have one or more _children_ containers. In the case of [Residue] instances, for example, consider the sequence _ALA-GLC-PRO_. In such a molecule, _ALA_ would be the _parent_ of _GLC_ which, in turn, would be the _parent_ of _PRO_. The same logic applies to relationships between atoms. This structural organization allows ProtoSyn to then infer parenthood up to _N_ levels. When conjugated with the internal coordinates system, this allows for easy and fast manipulation of dihedral angles, for example.

The information regarding these internal coordinates (and cartesian coordinates) is organized in the pose's **State**. This is subdivided into a list of [AtomState] instances and a [StateMatrix]. Both of these structures are interchangable, meaning that altering a value in [StateMatrix] updates the corresponding value in the correct [AtomState] and vice-versa. The rational behind having both structures lies in having the ease/speed of changing a large volume of coordinates in the [StateMatrix] at once, while still being able to control internal coordinates in single atoms using the [AtomState].

In order to initialize the internal coordinate system, each [Topology] has an extra set of 3 pseudo atoms, called a **Root**, that sets-up the ascendents for the first few atoms of the molecule.

## Drivers

![organization](../assets/ProtoSyn-drivers.png)

ProtoSyn.jl offers some quick simulation functionalities, such as [MonteCarlo], [ILS] and [SteepestDescent] simulations, using [Driver] instances. A Driver is a function which _drives_ the pose from one state to the next. As a general rule, these Drivers are usually comprised of two important components: [Mutator] instances and an [EnergyFunction].

A [Mutator] is a function that performs a single change in the system. For example, a [DihedralMutator] will rotate a random dihedral by a random amount. These can be parametrized for more specific needs.

An [EnergyFunction] evaluates the fitness of a given Pose based on a set of [EnergyFunctionComponent] instances. These can be distance-based restrictions, positional agreement with machine learning models (such as [TorchANI]), etc. 

!!! ukw "ProtoSyn tip"
    Some commonly used [EnergyFunction]s can be found in the [Common] module.

## Modular system

ProtoSyn.jl package is organized in a modular fashion, with each new module adding specific methods and features relative to a given scientific area of expertise. All these modules build on top of the base **Core** module. For example, in the Core module, [ProtoSyn.load] is able to read a PDB file into a [Pose] struct, but calling [Peptides.load] \(from the [Peptides] module) will **also** infer peptidic connections between amino acids, setting the correct parenthood relationships between residues in the pose's graph.

# Next steps

* We recommend you check the [Examples] page for some initial tutorials on how to use some of ProtoSyn's functionalities.
* If you want a deeper dive into the inner workings of this package, check the [ProtoSyn.jl API] page. 