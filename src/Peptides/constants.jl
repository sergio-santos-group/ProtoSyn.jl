using CSV

# Positive values: Hydrophobic
# Negative values: Hydrophilic
const doolitle_hydrophobicity = Dict{String, ProtoSyn.Units.defaultFloat}(
    "ILE" =>  4.5, # I
    "VAL" =>  4.2, # V
    "LEU" =>  3.8, # L
    "PHE" =>  2.8, # F
    "CYS" =>  2.5, # C
    "MET" =>  1.9, # M
    "ALA" =>  1.8, # A
    "GLY" => -0.4, # G
    "THR" => -0.7, # T
    "SER" => -0.8, # S
    "TRP" => -0.9, # W
    "TYR" => -1.3, # Y
    "PRO" => -1.6, # P
    "HIS" => -3.2, # H
    "HIE" => -3.2, # H
    "HID" => -3.2, # H
    "ASN" => -3.5, # N
    "GLU" => -3.5, # E
    "GLN" => -3.5, # Q
    "ASP" => -3.5, # D
    "LYS" => -3.9, # K
    "ARG" => -4.5  # R
)

const doolitle_hydrophobicity_mod1 = Dict{String, ProtoSyn.Units.defaultFloat}(
    "ILE" =>  5.5,
    "VAL" =>  5.2,
    "LEU" =>  4.8,
    "PHE" =>  3.8,
    "CYS" =>  3.5,
    "MET" =>  2.9,
    "ALA" =>  2.8,
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

const doolitle_hydrophobicity_mod2 = Dict{String, ProtoSyn.Units.defaultFloat}(
    "ILE" =>  6.5,
    "VAL" =>  6.2,
    "LEU" =>  5.8,
    "PHE" =>  4.8,
    "CYS" =>  4.5,
    "MET" =>  3.9,
    "ALA" =>  3.8,
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


const doolitle_hydrophobicity_extreme = Dict{String, ProtoSyn.Units.defaultFloat}(
    "ILE" => 11.5,
    "VAL" => 11.2,
    "LEU" => 10.8,
    "PHE" => 9.8,
    "CYS" => 9.5,
    "MET" => 8.9,
    "ALA" => 8.8,
    "GLY" => 0.0,
    "THR" => 0.0,
    "SER" => 0.0,
    "TRP" => 0.0,
    "TYR" => 0.0,
    "PRO" => 0.0,
    "HIS" => 0.0,
    "HIE" => 0.0,
    "ASN" => 0.0,
    "GLN" => 0.0,
    "ASP" => 0.0,
    "GLU" => 0.0,
    "LYS" => 0.0,
    "ARG" => 0.0
)

available_aminoacids = Dict{Char, Bool}()

const polar_residues = ["ARG", "ASN", "ASP", "GLU", "GLN", "HIS", "LYS", "SER", "THR"]

function load_aa_similarity(::Type{T}, filename::String) where {T <: AbstractFloat}

    _filename = joinpath(Peptides.resource_dir, filename)
    data = CSV.File(_filename; skipto = 2)

    aa_similarity = Dict{String, Dict{String, T}}()

    for row in data
        k = string.(keys(row)[2:end])
        aa_similarity[row[:X]] = Dict(zip(k, values(row)[2:end]))
    end

    return aa_similarity
end

const aminoacid_similarity = load_aa_similarity(ProtoSyn.Units.defaultFloat, "aa_similarity.csv")