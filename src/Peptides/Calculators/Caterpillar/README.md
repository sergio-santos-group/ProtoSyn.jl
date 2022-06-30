# Caterpillar solvation energy

The Caterpillar solvation energy calculation is based on the work of Coluzza
et al (see https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0020853).

The ProtoSyn's modifications introduce a significant degree of complexity, 
futher explained bellow. The calculation takes 3 steps:


## 1. Burial degree calculation

In this step, a given algorithm loops over all the residues (selecting a given
atom for distance matrix calculation) and identifies the burial degree of each
residue (this step can be parametrized by (1) the burial degree algorithm, (2)
the identification curve, (3) the selection atom, (4) the rmax cut-off and (5)
slope control sc. These settings are further explained bellow).

### 1.1 Burial degree algorithm
ProtoSyn offers 2 different burial degree identification algorithms, the
Neighbour Count (NC) and Neighbour Vector (NV) (as explained further in
https://pubmed.ncbi.nlm.nih.gov/19234730/). In comparison to eachother, NC
algorithms only take into consideration the number of selected atoms within a
defined rmax range, while NV algorithms also take into consideration the
orientation of neighbouring residues to defined the burial degree.

### 1.2 Identification curve
The available identification curves are linear, sigmoid and normalized sigmoid
(in NV algorithms only). Note that the definition of the identification curve
controls the amount of distance information considered for the calculation of
the burial degree: linear identification curves incorporate more distance
information than sigmoid identification curves (normalized sigmoid, in NV
algorithms, use the least ditance information, similar to NC algorithms).

### 1.3 Selection atom
This can be any selection that yields an atom. However, in most cases applied
to proteins, either the Cα or Cβ atoms should be chosen.

### 1.4 rmax cutoff
This can be any float number, however, the range between 9.0Å and 50.0Å was
identified as yielding the best results. Note that short rmax values cause a
much more localized identification of burial degrees (i.e.: in comparison with
the rest of the secondary structure) while larger rmax values identify a more
global level of burial (i.e.: in comparison with all aminoacids in the
structure).

### 1.5 slope control (sc)
This value controls the slope degree in sigmoid identification curves. Lower
values yield less pronounced slopes (therefore taking more distance
information into consideration) while higher sc values define more strict
cut off lines.