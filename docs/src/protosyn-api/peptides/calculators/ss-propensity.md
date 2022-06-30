```@meta
CurrentModule = ProtoSyn.Peptides.Calculators
```

# Secondary Structure Propensity

“Propensity for secondary structures represents an intrinsic property of an amino acid (…)” - Susan Contantini, 2006

The [Secondary Structure Propensity](@ref) Calculators module introduces a measure of how likely a given peptide sequence is based on the natural distribution of aminoacids in nature for the assessed secondary structure, according to Constantini et al. (See [this paper](https://www.sciencedirect.com/science/article/pii/S0006291X06002543)).

```@docs
load_default_aa_ss_propensity
calc_aa_ss_propensity
get_default_aa_ss_propensity
fixate_secondary_structure!
```