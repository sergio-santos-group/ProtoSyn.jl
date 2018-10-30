# Input JSON Schema

This section describes in detail the general schematics of the input JSON so it can be easily read by ProtoSyn. 

## Topology

Protosyn receives various types of topology, depending on the [Driver](@ref Drivers) chosen by the user.

### Forcefield Topology

Forcedield topology holds the information regarding the system bonds, angles, dihedrals and non-bonded interactions. ProtoSyn offers a tool ([ProtoSyn.jl/tools/tpr2json.py](https://github.com/sergio-santos-group/ProtoSyn.jl/blob/master/tools/tpr2json.py)) that easily reads .tpr file formats from [GROMACS](http://www.gromacs.org) (dumped to human readable format using `gmx dump`) and converts the obtained information to a correct JSON format. 

```@raw html
<details>
<summary style="background-color:#1abc9c;border-radius: 13px;height: 30px;display: flex;align-items: center">&ensp;Detailed forcefield topology JSON description▼</summary>
<p>
```
```javascript
{
    "dihedraltypes": {
        ID::String: {
            "phi"    : Float64,    // Equilibrium angle value (degrees)
            "cp"     : Float64,    // Force constant (kJ mol⁻¹)
            "mult"   : Float64     // Multiplicity
        },
        ...
    },
    "angletypes": {
        ID::String: {
            "th"     : Float64,    // Equilibrium angle value (degrees)
            "ct"     : Float64     // Force constant (kJ mol⁻¹ rad⁻²)
        },
        ...
    },
    "bondtypes": {
        ID::String: {
            "cb"     : Float64,    // Force constant (kJ mol⁻¹ nm⁻²)
            "b0"     : Float64     // Equilibrium bond length (nm)
        },
        ...
    },
    "atomtypes": {
        ID::String: {
            "sigma"  : Float64,    // Finite distance at which the inter-particle potential is zero (nm)
            "epsilon": Float64,    // Depth of the potential well (kJ mol⁻¹)
            "name"   : String      // Name of the atomtype
        },
        ...
    },
    "dihedrals": [
        {
            "a1"     : Int64,      // Global atom 1 index
            "a2"     : Int64,      // Global atom 2 index
            "a3"     : Int64,      // Global atom 3 index
            "a4"     : Int64,      // Global atom 4 index
            "type"   : Int64       // ID of dihedraltypes
        },
        ...
    ],
    "angles": [
        {
            "a1"     : Int64,      // Global atom 1 index
            "a2"     : Int64,      // Global atom 2 index
            "a3"     : Int64,      // Global atom 3 index
            "type"   : Int64       // ID of angletypes
        },
        ...
    ],
    "bonds": [
        {
            "a1"     : Int64,      // Global atom 1 index
            "a2"     : Int64,      // Global atom 2 index
            "type"   : Int64       // ID of bondtypes
        },
        ...
    ]
    "atoms": [
        {
            "q"      : Float64,    // Atom charge (eletron)
            "type"   : Int64,      // ID of atomtypes
        },
        ...
    ],
    "excls": {
        a1::Int64: Array{Int64, 1}, // Global indices of the atoms excluded from non-bonded interactions with a1
        ...
    },
    "pairs": [
        {
            "a1"     : Int64,       // Global index of a1
            "a2"     : Int64,       // Global index of a2 (At 3 connections away from a1)
        },
        ...
    ]
}
```
```@raw html
</p>
</details>
```

### Monte Carlo Topology

Monte Carlo topology holds the information regarding the rotatable dihedrals and the residues present in the system. This allows for correct rotation in Monte Carlo moves.

```@raw html
<details>
<summary style="background-color:#1abc9c;border-radius: 13px;height: 30px;display: flex;align-items: center">&ensp;Detailed Monte Carlo JSON description▼</summary>
<p>
```
```javascript
{
    "dihedrals": [
        {
            "parent" : Int64,           // Index of parent residue in list of residues
            "movable": Array{Int64, 1}, // Global indices of atoms moved in the parent residue when rotating this dihedral
            "a1"     : Int64,           // Global atom 1 index
            "a2"     : Int64,           // Global atom 2 index
            "a3"     : Int64,           // Global atom 3 index
            "a4"     : Int64,           // Global atom 4 index
            "type"   : String           // Dihedral type (i.e. "PHI", "PSI", ...)
        },
        ...
    ],
    "residues": [
        {
            "next"   : Int64,           // Index of the next residue in list of residues
            "atoms"  : Array{Int64, 1}, // List of global indices of the atoms included in this residue
            "type"   : String,          // Residue type (i.e. "V", "K", ...)
            "n"      : Int64            // Residue identifier
        },
        ...
    ]
}
```
```@raw html
</p>
</details>
```