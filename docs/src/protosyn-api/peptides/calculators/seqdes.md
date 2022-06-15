```@meta
CurrentModule = ProtoSyn.Peptides.Calculators.SeqDes
```

# SeqDes model

The [SeqDes model](@ref) Calcualtors module introduces the possibility to employ part of the work by Huang et al. (see [this paper](https://www.nature.com/articles/s41467-022-28313-9)) in ProtoSyn. The [SeqDes model](@ref), in sum, evaluates the current sequence and chi dihedrals in accordance to a trained model on many natural proteins (but generalizable to de novo sequences).

```@docs
get_pdb_data
calc_seqdes
get_default_seqdes
```