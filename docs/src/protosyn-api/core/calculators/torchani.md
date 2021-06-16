```@meta
CurrentModule = ProtoSyn.Calculators
```

# TorchANI

[TorchANI](https://github.com/aiqm/torchani) is a Python implementation of [ANI](https://github.com/isayev/ASE_ANI) machine learning model for energy calculation. Trained on a large dataset of DFT calculation results, TorchANI is able to obtain DFT-level accuracy in energy calculation in a fraction of the time. For more details, read the [original paper](https://pubs.acs.org/doi/10.1021/acs.jcim.0c00451). Making use of the [PyCall](https://github.com/JuliaPy/PyCall.jl) package and Julia's native methods for code integration, ProtoSyn incorporates TorchANI as an [`EnergyFunctionComponent`](@ref).

```@docs
Calculators.TorchANI.get_ani_species
Calculators.TorchANI.calc_torchani_ensemble
Calculators.TorchANI.calc_torchani_model
Calculators.TorchANI.get_default_torchani_ensemble
Calculators.TorchANI.get_default_torchani_model
```

![ProtoSyn TorchANI Components](../../../assets/ProtoSyn-torchani-comp.png)

**Figure 1 |** A diagram representation of the default `TorchANI ML Ensemble` and `TorchANI ML Model` [`EnergyFunctionComponent`](@ref) instances, obtained by usign the [`get_default_torchani_ensemble`](@ref Calculators.TorchANI.get_default_torchani_ensemble) and [`get_default_torchani_model`](@ref Calculators.TorchANI.get_default_torchani_model) methods, respectively. While the `TorchANI ML Ensemble` is more accurate, it's also significantly slower. Depending on the application, employing a single model of the ensemble (of 8 models) might be sufficient.

ProtoSyn also makes available the `TorchANI ML Model` [`EnergyFunctionComponent`](@ref) instance as a call to an [XML-RPC server](https://en.wikipedia.org/wiki/XML-RPC). This protocol spawns a Python server in parallel, running TorchANI, who receives XML requests from Julia and returns the calculated energy and forces. This functionality might be useful in certain systems and machines.

```@docs
Calculators.TorchANI.start_torchANI_server
Calculators.TorchANI.stop_torchANI_server
Calculators.TorchANI.calc_torchani_model_xmlrpc
Calculators.TorchANI.get_default_torchani_model_xmlrpc
```