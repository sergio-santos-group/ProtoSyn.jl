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

baremodule Dihedral
    struct DihedralType
        name::String
        atom::String
    end

    const psi   = DihedralType("psi", "N")
    const ψ     = DihedralType("psi", "N")

    const omega = DihedralType("omega", "CA")
    const ω     = DihedralType("omega", "CA")
    
    const phi   = DihedralType("phi", "C")
    const ϕ     = DihedralType("phi", "C")
end


export SecondaryStructure
# phi, psi, omega
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