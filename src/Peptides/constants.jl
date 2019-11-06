const one_2_three = Dict{Char,String}(
    '?' => "BKB",
    'A' => "ALA",
    'C' => "CYS",
    'D' => "ASP",
    'E' => "GLU",
    'F' => "PHE",
    'G' => "GLY",
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

const three_2_one = Dict{String,Char}(v=>k for (k,v) in one_2_three)


baremodule DihedralTypes
using Base: @enum
# using Base
@enum Type begin
    undef = 0
    phi   = 1
    psi   = 2
    omega = 3
    chi1  = 4
    chi2  = 5
    chi3  = 6
    chi4  = 7
    chi5  = 8
    # ϕ = 1
    # ψ = 2
    # ω = 3
    # χ1 = 4
    # χ2 = 5
    # χ3 = 6
    # χ4 = 7
    # χ5 = 8

end

# Base.show(io::IO, item::Type) = begin
#     Base.println(io, "iiiee")
# end

end



# struct SecondaryStructureSpec
#     name::String
#     atoms::Tuple{Vararg{String}}
# end

# const ss2bbd = Dict(
#         Common.SS.SHEET => Dict(Common.DIHEDRAL.phi => deg2rad(-135.0), Common.DIHEDRAL.psi => deg2rad(135.0)),
#         Common.SS.HELIX => Dict(Common.DIHEDRAL.phi => deg2rad(-60.0),  Common.DIHEDRAL.psi => deg2rad(-45.0)))
# end

# const SecondaryStructure = (
#     helix=()
# )