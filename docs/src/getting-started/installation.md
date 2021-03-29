# Installation

ProtoSyn.jl version 1.0 has been developed in tested for [Julia 1.5](https://julialang.org/downloads/).

After setting up you Julia installation, open a new REPL and add ProtoSyn.jl using the package manager.

```@julia
julia> ] add ProtoSyn.jl
```

All Julia-based dependencies should automatically installed/updated. In order to use TorchANI's energy function, you'll also need to install [Python](https://www.python.org/downloads/) with both [Torch](https://pytorch.org/get-started/locally/) and [TorchANI](https://aiqm.github.io/torchani/start.html) libraries.

A visualization macro for [PyMOL](https://pymol.org/2/) is included in ProtoSyn's core module, and it's the reccomended visualization tool. 