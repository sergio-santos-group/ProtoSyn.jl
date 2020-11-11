const one_2_three = Dict{Char,String}(
    '?' => "BKB",
    'A' => "ALA",
    'C' => "CYS",
    'D' => "ASP",
    'E' => "GLU",
    'F' => "PHE",
    'G' => "GLY",
    'H' => "HIS",
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

const three_2_one = Dict{String, Char}(v=>k for (k,v) in one_2_three)


module Dihedral

    using ProtoSyn: Residue

    abstract type DihedralType end

    # --- Backbone dihedral angles

    struct Phi <: DihedralType end; const phi = Phi()
    struct Psi <: DihedralType end; const psi = Psi()
    struct Omega <: DihedralType end; const omega = Omega()

    (phi::Phi)(residue::Residue)     = residue["C"]
    (psi::Psi)(residue::Residue)     = residue["N"]
    (omega::Omega)(residue::Residue) = residue["CA"]

    # --- Chi dihedral angles

    struct Chi1 <: DihedralType end; const chi1 = Chi1()
    struct Chi2 <: DihedralType end; const chi2 = Chi2()
    struct Chi3 <: DihedralType end; const chi3 = Chi3()
    struct Chi4 <: DihedralType end; const chi4 = Chi4()

    function get_chi(residue::Residue, chi::Int)
        available_chis = chi_dict[residue.name]
        @assert chi <= length(available_chis) "Tried to retrieve chi $chi on a residue with only $(length(available_chis)) chi angles defined"
        chi_atom = available_chis[chi]
        @assert residue[chi_atom] !== nothing "Chi $chi of residue $residue requires atom $(chi_atom.name), which was not found"
        @assert length(residue[chi_atom].children) > 0 "At least one children atom of $(residue[chi_atom]) needs to exist (0 found)."
        
        return residue[chi_atom].children[1]
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
    const chi_dict = Dict{String, Vector{String}}(
        "ALA" => [],
        "CYS" => ["CB"],
        "ASP" => ["CB", "CG"],
        "GLU" => ["CB", "CG", "CD"],
        "PHE" => ["CB", "CG"],
        "GLY" => [],
        "HIS" => ["CB", "CG"],
        "ILE" => ["CB", "CG1"],
        "LYS" => ["CB", "CG", "CD", "CE"],
        "LEU" => ["CB", "CG"],
        "MET" => ["CB", "CG", "SD"],
        "ASN" => ["CB", "CG"],
        "PRO" => ["CB", "CG", "CD"], # Attention
        "GLN" => ["CB", "CG", "CD"],
        "ARG" => ["CB", "CG", "CD", "NE"],
        "SER" => ["CB"],
        "THR" => ["CB"],
        "VAL" => ["CB"],
        "TRP" => ["CB", "CG"],
        "TYR" => ["CB", "CG"]
    )
end


export SecondaryStructure
#                                           phi,   psi, omega
const SecondaryStructure = Dict{Symbol, NTuple{3, Number}}(
    :antiparallel_sheet => map(deg2rad, (-139.0, 135.0, 180.0)),
    :parallel_sheet     => map(deg2rad, (-119.0, 113.0, 180.0)),
    :linear             => map(deg2rad, ( 180.0, 180.0, 180.0)),
    :helix              => map(deg2rad, ( -60.0, -45.0, 180.0))
)

const doolitle_hydrophobicity = Dict{String, ProtoSyn.Units.defaultFloat}(
    "ALA" =>  1.8,
    "CYS" =>  2.5,
    "ASP" => -3.5,
    "GLU" => -3.5,
    "PHE" =>  2.8,
    "GLY" => -0.4,
    "HIS" => -3.2,
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

const doolitle_hydrophobicity_mod = Dict{String, ProtoSyn.Units.defaultFloat}(
    "ILE" =>  10.5,
    "VAL" =>  10.2,
    "LEU" =>  9.8,
    "PHE" =>  8.8,
    "CYS" =>  8.5,
    "MET" =>  7.9,
    "ALA" =>  7.8,
    "GLY" => -0.4,
    "THR" => -0.7,
    "SER" => -0.8,
    "TRP" => -0.9,
    "TYR" => -1.3,
    "PRO" => -1.6,
    "HIS" => -3.2,
    "ASN" => -3.5,
    "GLN" => -3.5,
    "ASP" => -3.5,
    "GLU" => -3.5,
    "LYS" => -3.9,
    "ARG" => -4.5
)