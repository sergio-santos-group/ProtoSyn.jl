```@meta
CurrentModule = ProtoSyn.Calculators.GB
```

# Generalized Born Solvation

The Generalized Born (GB) solvation model attempts to approximate the Poisson-Boltzmann equation, using the following equation:

$$- \frac{1}{2}\left ( \frac{1}{\varepsilon_{protein} }-\frac{1}{\varepsilon_{solvent} } \right ) \sum_{i,j}^{N}\frac{q_{i}q_{j}}{f_{GB}},$$

$$f_{GB} = \sqrt{d_{ij}^{2} + \left ( \alpha_{i}\alpha_{j} e^{-\frac{d_{ij}}{4\alpha_{i}\alpha_{j}}} \right)}$$

where $\varepsilon_{protein}$ and $\varepsilon_{solvent}$ are the dieletric constants of the protein interior and of the solvent, respectively, $qi$ and $qj$ are the atomic partial charges for each interacting point-like particle, $d_{ij}$ is the inter-atomic distance of the considered particle pair and, finally, $\alpha_{i}$ and $\alpha_{j}$ are the _effective Born radii_ of both considered interacting particles. The Born radius of an [`Atom`](@ref) is a characterization of its burial degree in the solute, and defines the interaction amount with the solvent. Despite the GB solvation model success, accurate estimation of the Born radii in a system has been the bottleneck preventing mass adoption as the go-to implicit solvent potential.

The following section has been subdivided in several sub-section for organizational purposes:

+ [Born Radii estimation](@ref)
+ [Generalized Born solvation EnergyFunctionComponent](@ref)

## Born Radii estimation

As previously stated, Born radii estimation is a performance-intensive task. Several approaches have been devised over the years to improve both the accuracy and calculation speed in Born radii estimation. ProtoSyn employs the work of Fogolari et al. (2020) (https://pubmed.ncbi.nlm.nih.gov/31693089/), estimating the Born radii using a machine learning model.

```@docs
predict_igbr_nn_born_radii
```

## Generalized Born solvation EnergyFunctionComponent

With the estimated Born radii, ProtoSyn employs the above-defined equations to evaluate a system's implicit solvation status using the following methods:

```@docs
calc_gb
get_default_gb
```

!!! ukw "Note:"
    Despite using modern machine learning models (which substantially cut down the performance cost of estimating accurate born radii), this step is still costly and time-consuming. In not big changes are expected in the system (for example, only small refinement movements are introduced in the system, such as when re-packaging sidechains), consider using static born radii (calculated once and provided to [`calc_gb`](@ref) as a static `Vector`).

!!! ukw "Note:"
    Modern implicit solvation models often employ a hybrid approach, dubbed "SASA/GB". In short, the Generalized Born model attempts to estimate the enthalpic contribution of solvating a given molecule, and the SASA model calculates the entropy contribution of "opening space" for the solvation of such a molecule. Both models act together to provide a more clear picture of the solvation potential of that system. Consider employing the [`get_default_sasa`](@ref ProtoSyn.Calculators.SASA.get_default_sasa) [`EnergyFunctionComponent`](@ref) in conjunction with the [`get_default_gb`](@ref ProtoSyn.Calculators.GB.get_default_gb) [`EnergyFunctionComponent`](@ref).