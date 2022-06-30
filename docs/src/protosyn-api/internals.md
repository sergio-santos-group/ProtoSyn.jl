# Internals

## How did you get here?!

The following docstrings are a dump for non-exported private & internal functions. These provide basic internal functionality that should rarely be employed by an end-user. In either case, docstrings exist to explain the inner-workings of methods for potential developers.

```@meta
CurrentModule = ProtoSyn
```

```@docs
ProtoSyn.Calculators.TorchANI.r_xml_travel!
ProtoSyn.Clustering.dunn_index
ProtoSyn.Calculators.get_available_energy_function_components
Base.resize!(::ProtoSyn.Calculators.VerletList, ::Int64)
ProtoSyn.Common
ProtoSyn.AbstractContainer
ProtoSyn.Peptides.Calculators.Caterpillar.nc_scalling_exposed_only
ProtoSyn.Peptides.Calculators.Caterpillar.nv_scalling_exposed_only
ProtoSyn.Peptides.Calculators.Caterpillar.nc_non_scalling_exposed_only
ProtoSyn.Peptides.Calculators.Caterpillar.nv_non_scalling_exposed_only
ProtoSyn.Peptides.Calculators.Caterpillar.nc_scalling_all_contributions
ProtoSyn.Peptides.Calculators.Caterpillar.nv_scalling_all_contributions
ProtoSyn.Peptides.Calculators.Caterpillar.nc_non_scalling_all_contributions
ProtoSyn.Peptides.Calculators.Caterpillar.nv_non_scalling_all_contributions
ProtoSyn.cross2d
ProtoSyn.@cross
ProtoSyn.@dot
ProtoSyn.Calculators.@reduce
ProtoSyn.get_lead
ProtoSyn.Calculators.distance_matrix_kernel
ProtoSyn.Clustering.complete_diameter_distances
ProtoSyn.Clustering.davies_bouldin_index
ProtoSyn.Clustering.rmsd_matrix
ProtoSyn.Clustering.average_diameter_distances
ProtoSyn.tile
ProtoSyn.tile!
ProtoSyn.Mutators
ProtoSyn.Peptides.Mutators
ProtoSyn.Peptides
ProtoSyn.Materials
ProtoSyn.GMX
ProtoSyn.Peptides.GMX
ProtoSyn.GMX.check_installation
ProtoSyn.Calculators.neighbours
ProtoSyn.Calculators.Restraints
ProtoSyn.print_diagnose_results
ProtoSyn.XMLRPC.ClientProxy
ProtoSyn.read_yml
ProtoSyn.gpu_allocation
ProtoSyn.promote
```