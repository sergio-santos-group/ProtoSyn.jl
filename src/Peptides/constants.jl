const one_2_three = Dict{Char,String}(
    '?' => "BKB",
    'A' => "ALA",
    'C' => "CYS",
    'D' => "ASP",
    'E' => "GLU",
    'F' => "PHE",
    'G' => "GLY",
    'H' => "HIS",
    'H' => "HIE",
    'I' => "ILE",
    'K' => "LYS",
    'L' => "LEU",
    'M' => "MET",
    'N' => "ASN",
    'P' => "PRO",
    'Q' => "GLN",
    'R' => "ARG",
    'S' => "SER",
    'T' => "THR",
    'V' => "VAL",
    'W' => "TRP",
    'Y' => "TYR",
)

const three_2_one = Dict{String, Char}(
    "BKB" => '?',
    "ALA" => 'A',
    "CYS" => 'C',
    "ASP" => 'D',
    "GLU" => 'E',
    "PHE" => 'F',
    "GLY" => 'G',
    "HIS" => 'H',
    "HIE" => 'H',
    "ILE" => 'I',
    "LYS" => 'K',
    "LEU" => 'L',
    "MET" => 'M',
    "ASN" => 'N',
    "PRO" => 'P',
    "GLN" => 'Q',
    "ARG" => 'R',
    "SER" => 'S',
    "THR" => 'T',
    "VAL" => 'V',
    "TRP" => 'W',
    "TYR" => 'Y',
)


"""
    Dihedral

An auxiliary struct for [`Dihedral`](@ref) selection by name. Each of the
available dihedrals (under the abstract type `DihedralType`, see bellow) can be
called using the following signature, returning the respective representative
atom:

```julia
(dihedral::DihedralType)(residue::Residue)
```

# Available dihedrals

`Dihedral.phi` `Dihedral.psi` `Dihedral.omega`
`Dihedral.chi1` `Dihedral.chi2` `Dihedral.chi3` `Dihedral.chi4`

# See also
[`getdihedral`](@ref ProtoSyn.getdihedral)
[`setdihedral!`](@ref ProtoSyn.setdihedral!)

# Examples
```jldoctest
julia> ProtoSyn.Peptides.Dihedral.phi(pose.graph[1][1])
Atom{/UNK:1/UNK:1/SER:1/C:10}

julia> atom = ProtoSyn.Peptides.Dihedral.chi2(pose.graph[1][2])
Atom{/UNK:1/UNK:1/GLU:2/CD:22}

julia> ProtoSyn.getdihedral(pose.state, atom)
-3.141165
```

"""
module Dihedral

    using ProtoSyn: Residue

    abstract type DihedralType end

    # --- Backbone dihedral angles

    struct Phi <: DihedralType end; const phi = Phi()
    struct Psi <: DihedralType end; const psi = Psi()
    struct Omega <: DihedralType end; const omega = Omega()

    (phi::Phi)(residue::Residue)     = residue["C"]
    (psi::Psi)(residue::Residue)     = begin
        if length(residue.children) == 0 
            @warn "Residue $residue has no psi angle"
            return nothing
        end
        residue.children[1]["N"]
    end
    (omega::Omega)(residue::Residue) = residue["CA"]

    # --- Chi dihedral angles

    struct Chi1 <: DihedralType end; const chi1 = Chi1()
    struct Chi2 <: DihedralType end; const chi2 = Chi2()
    struct Chi3 <: DihedralType end; const chi3 = Chi3()
    struct Chi4 <: DihedralType end; const chi4 = Chi4()

    function get_chi(residue::Residue, chi::Int)
        available_chis = chi_dict[residue.name.content]
        @assert chi <= length(available_chis) "Tried to retrieve chi $chi on a residue with only $(length(available_chis)) chi angles defined"
        chi_atom = available_chis[chi + 1] # Note the "+ 1"
        @assert residue[chi_atom] !== nothing "Chi $chi of residue $residue requires atom $(chi_atom), which was not found"
        
        return residue[chi_atom]
    end

    (chi1::Chi1)(residue::Residue) = get_chi(residue, 1)
    (chi2::Chi2)(residue::Residue) = get_chi(residue, 2)
    (chi3::Chi3)(residue::Residue) = get_chi(residue, 3)
    (chi4::Chi4)(residue::Residue) = get_chi(residue, 4)

    # Dict built according to Dunbrack rotamer library.
    # Dihedral angles not in the dict can still be rotated, by manually defining
    # the controling atom. For example, ALA-CB dihedral can be rotated by
    # setting the dihedral value Δϕ on any of its children. By definition, 
    # dihedral angles which only rotate the position of hydrogens are not
    # considered in this dict.
    chi_dict = Dict{String, Vector{String}}(
        "ALA" => ["HB1"],
        "CYS" => ["CB", "SG"],
        "ASP" => ["CB", "CG", "OD1"],
        "GLU" => ["CB", "CG", "CD", "OE1"],
        "PHE" => ["CB", "CG", "CD1"],
        "GLY" => [],
        "HIS" => ["CB", "CG", "ND1"],
        "HIE" => ["CB", "CG", "ND1"],
        "ILE" => ["CB", "CG1", "CD", "HD1"],
        "LYS" => ["CB", "CG", "CD", "CE", "NZ"],
        "LEU" => ["CB", "CG", "CD1"],
        "MET" => ["CB", "CG", "SD", "CE"],
        "ASN" => ["CB", "CG", "OD1", "HD21"],
        "PRO" => ["CB", "CG", "CD", "HD1"], # Attention
        "GLN" => ["CB", "CG", "CD", "NE2"],
        "ARG" => ["CB", "CG", "CD", "NE", "CZ"],
        "SER" => ["CB", "OG"],
        "THR" => ["CB", "OG1"],
        "VAL" => ["CB", "CG1"],
        "TRP" => ["CB", "CG", "CD1"],
        "TYR" => ["CB", "CG", "CD1"]
    )
end

# --- Secondary Structure ------------------------------------------------------

export SecondaryStructure
export SecondaryStructureTemplate

"""
    WIP
"""
struct SecondaryStructureTemplate{T <: AbstractFloat}
    ϕ::T
    ψ::T
    ω::T
end

SecondaryStructureTemplate(template::Tuple{T, T, T}) where {T <: AbstractFloat} = begin
    return SecondaryStructureTemplate(template[1], template[2], template[3])
end

"""
    WIP
"""
const SecondaryStructure = Dict{Symbol, SecondaryStructureTemplate}(
    #                                                                  phi,   psi, omega
    :antiparallel_sheet => SecondaryStructureTemplate(map(deg2rad, (-130.0, 145.0, 180.0))),
    :parallel_sheet     => SecondaryStructureTemplate(map(deg2rad, (-110.0, 120.0, 180.0))),
    :linear             => SecondaryStructureTemplate(map(deg2rad, ( 180.0, 180.0, 180.0))),
    :helix              => SecondaryStructureTemplate(map(deg2rad, ( -60.0, -45.0, 180.0)))
)

const doolitle_hydrophobicity = Dict{String, ProtoSyn.Units.defaultFloat}(
    "ALA" =>  1.8,
    "CYS" =>  2.5,
    "ASP" => -3.5,
    "GLU" => -3.5,
    "PHE" =>  2.8,
    "GLY" => -0.4,
    "HIS" => -3.2,
    "HIE" => -3.2,
    "ILE" =>  4.5,
    "LYS" => -3.9,
    "LEU" =>  3.8,
    "MET" =>  1.9,
    "ASN" => -3.5,
    "PRO" => -1.6,
    "GLN" => -3.5,
    "ARG" => -4.5,
    "SER" => -0.8,
    "THR" => -0.7,
    "VAL" =>  4.2,
    "TRP" => -0.9,
    "TYR" => -1.3
)

const doolitle_hydrophobicity_mod3 = Dict{String, ProtoSyn.Units.defaultFloat}(
    "ILE" =>  7.5,
    "VAL" =>  7.2,
    "LEU" =>  6.8,
    "PHE" =>  5.8,
    "CYS" =>  5.5,
    "MET" =>  4.9,
    "ALA" =>  4.8,
    "GLY" => -0.4,
    "THR" => -0.7,
    "SER" => -0.8,
    "TRP" => -0.9,
    "TYR" => -1.3,
    "PRO" => -1.6,
    "HIS" => -3.2,
    "HIE" => -3.2,
    "ASN" => -3.5,
    "GLN" => -3.5,
    "ASP" => -3.5,
    "GLU" => -3.5,
    "LYS" => -3.9,
    "ARG" => -4.5
)

const doolitle_hydrophobicity_mod7 = Dict{String, ProtoSyn.Units.defaultFloat}(
    "ILE" =>  11.5,
    "VAL" =>  11.2,
    "LEU" =>  10.8,
    "PHE" =>  9.8,
    "CYS" =>  9.5,
    "MET" =>  8.9,
    "ALA" =>  8.8,
    "GLY" => -0.4,
    "THR" => -0.7,
    "SER" => -0.8,
    "TRP" => -0.9,
    "TYR" => -1.3,
    "PRO" => -1.6,
    "HIS" => -3.2,
    "HIE" => -3.2,
    "ASN" => -3.5,
    "GLN" => -3.5,
    "ASP" => -3.5,
    "GLU" => -3.5,
    "LYS" => -3.9,
    "ARG" => -4.5
)

available_aminoacids = ['M', 'K', 'P', 'Q', 'I', 'H', 'E', 'W', 'S', 'T', 'C', 'D', 'A', 'L', 'Y', 'V', 'R', 'G', 'F', 'N']