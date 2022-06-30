# Installation

ProtoSyn.jl version 1.1 has been developed and tested for [Julia 1.7.2](https://julialang.org/downloads/).

After setting up you Julia installation, open a new REPL and add ProtoSyn.jl using the package manager, directly from the project's GitHub page.

```@julia
julia> ] add https://github.com/sergio-santos-group/ProtoSyn.jl.git
```

All Julia-based dependencies should automatically be installed/updated.

In order to use some energy function components, such as TorchANI or REF15 (from PyRosetta), a [Python](https://www.python.org/downloads/) installation is required. ProtoSyn attempts to locate the necessary packages in the system (from `PATH` and `PYTHON_PATH`). If not available, installation instructions are shown for each missing package.

For trying out new features, consider using the development branch of ProtoSyn. Be aware that some bugs and missing documentation are to be expected when using the latest versions.

```@julia
julia> ] add https://github.com/sergio-santos-group/ProtoSyn.jl.git#dev
```