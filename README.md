<p align="center"> 
  <img src="./docs/src/assets/logo-white.png" alt="Logo">
</p>
<h3 align="center"> Sérgio M. Santos & José M. Pereira </h3>
<h5 align="center"> CICECO & Departamento de Química - <a href="https://www.ua.pt">Universidade de Aveiro</a></h5>


---


[![][docs-stable-img]][docs-stable-url] [![][build-status-img]][build-status-url]

## :scroll: Description

_ProtoSyn.jl is a platform for molecular manipulation and simulation, with an emphasis on peptide design._

The main goal of ProtoSyn is to be a basis on top of which new tools and protocols can be experimented and prototyped. Taking advantage of Julia’s environment, ProtoSyn has been built with emergent technologies in mind, such as distributed computing, GPU and SIMD acceleration and machine learning models usage.


## :round_pushpin: Features

* Edit a peptide structure by removing, adding and mutating any number of residues
* Create peptides from scratch by providing the desired sequence
* Copy parts of other molecules and graft them together to create something new
* Run Monte Carlo simulations to optimize a structure or a sequence
* Calculate energies and forces using [TorchANI](https://github.com/aiqm/torchani) code, integrated in Julia.
* Perform Steepest Descent optimizations
* Explore rotamer libraries to optimize sidechain packaging
* Select residues based on name, index, distance and other parameters, with a rich combinatory selection syntax
* Perform rigid body docking of ligands
* Include ramified carbohydrates and glycoproteins in your simulations, with support for sugar residues 


## :clipboard: Installation

After setting up you Julia installation, open a new REPL and add ProtoSyn.jl using the package manager.

```@julia
julia> ] add ProtoSyn.jl
```

All Julia-based dependencies should automatically installed/updated. In order to use TorchANI's energy function, you'll also need to install [Python](https://www.python.org/downloads/) with both [Torch](https://pytorch.org/get-started/locally/) and [TorchANI](https://aiqm.github.io/torchani/start.html) libraries.


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
