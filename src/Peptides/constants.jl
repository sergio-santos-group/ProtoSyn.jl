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
