<p align="center"> 
  <img src="./docs/src/assets/logo-white.png" alt="Logo">
</p>
<h3 align="center"> José M. Pereira & Sérgio M. Santos </h3>
<h5 align="center"> CICECO & Departamento de Química - <a href="https://www.ua.pt">Universidade de Aveiro</a></h5>

---

[![][docs-stable-img]][docs-stable-url] [![][build-status-img]][build-status-url]

## :scroll: Description

_ProtoSyn.jl is a platform for molecular manipulation and simulation, with an emphasis on peptide design._

The main goal of ProtoSyn.jl is to be a basis on top of which new tools and protocols can be experimented and prototyped. Taking advantage of Julia’s environment, ProtoSyn has been built with emergent technologies in mind, such as distributed computing, GPU/SIMD acceleration and machine learning models usage. ProtoSyn.jl intends to be a playground for trying out new emergent models and algorithms in a *plug-and-play* environment, completly open-source and with curated & up-do-date documentation, examples and tutorials.


## :round_pushpin: Features

* Edit a peptide structure by removing, adding and mutating any number of residues
* Create peptides from scratch by providing the desired sequence
* Copy parts of other molecules and graft them together to create something new
* Explore Ramachandran maps to introduce secondary structure variations 
* Run Monte Carlo simulations to optimize a structure or a sequence
* Calculate energies and forces using native potentials (such as SASA/GB, Coulomb or simple harmonics) or commonly used tools now integrated in Julia (such as [TorchANI](https://github.com/aiqm/torchani) or [REF15 from PyRosetta](https://www.rosettacommons.org/docs/latest/rosetta_basics/scoring/score-types)).
* Perform Steepest Descent optimizations
* Explore rotamer libraries to optimize sidechain packaging
* Select residues based on name, index, distance and other parameters, with a rich combinatory selection syntax
* Perform rigid body docking of ligands
* Include ramified carbohydrates and glycoproteins in your simulations, with support for sugar residues
* Include non-canonical aminoacids (NCAAs) and post-translational modifications, such as methylation or phosporylation
* Generate functionalized carbon models, including multi-layer & pore generation support


## :clipboard: Installation

After setting up your Julia installation, open a new REPL and add ProtoSyn.jl using the package manager:

```@julia
julia> ] add https://github.com/sergio-santos-group/ProtoSyn.jl.git
```

All Julia-based dependencies should be automatically installed/updated. Some functionalities of ProtoSyn.jl interact with external packages or tools (such as Python packages. For example, in order to use TorchANI's energy function, you'll also need to install [Python](https://www.python.org/downloads/) with both [Torch](https://pytorch.org/get-started/locally/) and [TorchANI](https://aiqm.github.io/torchani/start.html) libraries). In either case, upon using ProtoSyn.jl for the first time, users should be prompted to install any additional tool, if they so choose. These warning can be supressed at any time by adding `export JULIA_PROTOSYN_WARN_NON_AVALIABLE_EFC=false` in the user's `.bashrc` file.

For trying out new features, consider using the development branch of ProtoSyn.jl instead. Be aware that some bugs and missing documentation are to be expected when using the latest versions.

```@julia
julia> ] add https://github.com/sergio-santos-group/ProtoSyn.jl.git#dev
```

## :book: Documentation

- [**STABLE**][docs-stable-url] &mdash; documentation of the most recently tagged version.


## :email: Contacts

For any question or curiosity, please contact jose.manuel.pereira@ua.pt.

Check for the most recent developments in [ProtoSyn's blog.](https://sites.google.com/view/protosyn-jl/about)


## :trophy: Acknowledgments

<p align="center"> 
  <img src="./docs/src/assets/ProtoSyn-acknowledgments.png" alt="Logo">
</p>


[docs-stable-img]: https://img.shields.io/badge/docs-stable-blue.svg
[docs-stable-url]: https://sergio-santos-group.github.io/ProtoSyn.jl/stable

[build-status-img]: https://travis-ci.org/sergio-santos-group/ProtoSyn.jl.svg?branch=master
[build-status-url]: https://travis-ci.org/sergio-santos-group/ProtoSyn.jl
