@doc raw"""
    read_JSON(i_file::String)

Read a JSON file and return a dictionary with the content.

# Examples
```julia-repl
julia> Aux.read_JSON("i_file.json")
Dict{Any, Any}()
```
"""
function read_JSON(i_file::String)::Dict

    open(i_file, "r") do f
        json_txt = read(f, String)
        return JSON.parse(json_txt)
    end
end


@doc raw"""
    conv_21(aa::String)

Converts an aminoacid representation from 3 letters to 1, according to convention.

# Examples
```julia-repl
julia> Aux.conv321("VAL")
V
```
"""
function conv321(aa::String)

    conv321 = Dict(
    "CYS" => "C", "ASP" => "D", "SER" => "S", "GLN" => "Q", "LYS" => "K", "ILE" => "I", "PRO" => "P", "THR" => "T", "PHE" => "F", "ASN" => "N",
    "GLY" => "G", "HIS" => "H", "LEU" => "L", "ARG" => "R", "TRP" => "W", "ALA" => "A", "VAL" => "V", "GLU" => "E", "TYR" => "Y", "MET" => "M")
    if uppercase(aa) in keys(conv321)
        return conv321[uppercase(aa)]
    else
        return uppercase(aa)
    end
end


@doc raw"""
    conv123(aa::String)

Converts an aminoacid representation from 1 letter to 3, according to convention.

# Examples
```julia-repl
julia> Aux.conv321("V")
VAL
```
"""
function conv123(aa::String)

    conv123 = Dict(
    "Q" => "GLN", "W" => "TRP", "T" => "THR", "P" => "PRO", "C" => "CYS", "V" => "VAL", "L" => "LEU", "M" => "MET", "N" => "ASN", "H" => "HIS",
    "A" => "ALA", "D" => "ASP", "G" => "GLY", "E" => "GLU", "Y" => "TYR", "S" => "SER", "I" => "ILE", "K" => "LYS", "R" => "ARG", "F" => "PHE")
    if uppercase(aa) in keys(conv123)
        return conv123[uppercase(aa)]
    else
        return uppercase(aa)
    end
end