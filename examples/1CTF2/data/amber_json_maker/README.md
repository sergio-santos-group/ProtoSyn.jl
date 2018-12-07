# How to obtain a JSON file containing the Amber parameters from a PDB file?

## 1. Get a random MDP file and the pretended PDB file
This is a generic .mdp file for GROMACS. The contents of the file are irrelevant as long as it is a valid .mdp file.

```
random.mdp
1ctf.pdb
```

## 2. Place the pretended PDB file in a box
```bash
gmx editconf -f 1ctf.pdb -o 1ctf_boxed.pdb -c -box 4 -bt cubic
```

## 3. Get a valid topology
```bash
gmx pdb2gmx -f 1ctf.pdb -missing
python clean_topol.py > out.top
```

## 4. Create TPR
```bash
gmx grompp -f random.mdp -c 1ctf_boxed.pdb -p out.pdb -o 1ctf.tpr
```

## 5. Dump the TPR to a TXT file
```bash
gmx dump -s 1ctf.tpr > 1ctf.txt
```

## 6. Create JSON
```bash
python tpr2json.py 1ctf.txt > 1ctf_amber_top.json
```