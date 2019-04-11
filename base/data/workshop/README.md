```
     Created   José Pereira     March      2019
Last updated   José Pereira     March      2019
```

# Welcome to the workshop
Available tools:
- pdb_manipulator.py > Allows the manipulation of PDB files in order to regroup atom groups, reorder numbering in files, etc;
- clean_topol.py     > Cleans the GROMACS default topology file, searching for known irregularities;
- tpr2json.py        > Parses and outpus a readable JSON format from a dumped .txt TPR file;
- native2straight.jl > Reads a native structural PDB file and outputs the same structure with all PHI and PSI dihedral angles set to [-180, 180] degrees.

# Necessary files for ProtoSyn run

## 1. Amber JSON file
Refer to 'How to create ProtoSyn Input (2019).pdf' file. 

## 2. Straight protein structure PDF file
Refer to 'How to create ProtoSyn Input (2019).pdf' file. 

## 3. Secondary structure prediction
Refer to SCRATCH online server.

## 4. Contacts prediction
Refer to RaptorX online server.
