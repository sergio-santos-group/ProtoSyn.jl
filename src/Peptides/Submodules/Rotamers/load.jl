using ProtoSyn.Peptides: Dihedral

"""
    load_dunbrack([::Type{T}], [filename::String]) where {T <: AbstractFloat}

Load a Dunbrack styled rotamer library from file `filename` (if no `filename` is
given, will load the default rotamer library from the resources directory). If
no type `T` is provided, will use `ProtoSyn.Units.defaultFloat`. Return a
dictionary where each key is the name of an aminoacid type and the value is the
corresponding loaded [`BBD_RotamerLibrary`](@ref) instance. Note that not all
aminoacid types have sidechains, and therefore, have an associated 
[`BBD_RotamerLibrary`](@ref) instance.

# Examples
```
julia> rot_lib = ProtoSyn.Peptides.load_dunbrack(Float64)
Dict{String, ProtoSyn.Peptides.BBD_RotamerLibrary} with 19 entries:
  "GLN" => Name: GLN | Shape: (37, 37)…
  "LYS" => Name: LYS | Shape: (37, 37)…
  "ASN" => Name: ASN | Shape: (37, 37)…
  "TRP" => Name: TRP | Shape: (37, 37)…
  "THR" => Name: THR | Shape: (37, 37)…
  "VAL" => Name: VAL | Shape: (37, 37)…
  "HIS" => Name: HIS | Shape: (37, 37)…
  "SER" => Name: SER | Shape: (37, 37)…
  "PRO" => Name: PRO | Shape: (37, 37)…
  "ASP" => Name: ASP | Shape: (37, 37)…
  "PHE" => Name: PHE | Shape: (37, 37)…
  "ILE" => Name: ILE | Shape: (37, 37)…
  "TYR" => Name: TYR | Shape: (37, 37)…
  "HIE" => Name: HIS | Shape: (37, 37)…
  "ARG" => Name: ARG | Shape: (37, 37)…
  "LEU" => Name: LEU | Shape: (37, 37)…
  "MET" => Name: MET | Shape: (37, 37)…
  ⋮     => ⋮

```
"""
function load_dunbrack(::Type{T}, filename::String) where {T <: AbstractFloat}

    matrices = create_residue_library_matrices(T, filename)
    return fill_residue_library_matrices(T, filename, matrices)
end

function fill_residue_library_matrices(::Type{T}, filename::String, matrices::Dict{String, BBD_RotamerLibrary}) where {T <: AbstractFloat}

    open(filename) do file
        for line in readlines(file)
            startswith(line, "#") && continue
            elem = split(line)

            # Create a Rotamer from the current line
            #  a) Gather the name of the Rotamer corresponding residue
            name = string(elem[1])

            #  b) Gather chi values/standard deviation
            chis = Dict{DihedralType, Tuple{T, T}}()
            for index in 1:4
                elem[index + 4] == "0" && continue
                value = deg2rad(parse(T, elem[index + 9]))
                sd    = deg2rad(parse(T, elem[index + 13]))
                chis[getfield(Dihedral, Symbol("chi$index"))] = (value, sd)
            end

            #  c) Gather probability and create the Rotamer
            weight  = parse(T, elem[9])
            rotamer = Rotamer(name, chis)

            # Gather the location to append this rotamer to, in the rotamers
            # matrix
            rl = matrices[name]
            phi_index = findnearest(rl.phis, deg2rad(parse(T, elem[2])))
            psi_index = findnearest(rl.psis, deg2rad(parse(T, elem[3])))
        
            # Append to the found location. If necessary, initiate a new vector
            try
                push!(rl.rotamer_stacks[phi_index, psi_index], rotamer, weight)
            catch UndefRefError
                rl.rotamer_stacks[phi_index, psi_index] = RotamerStack(T)
                push!(rl.rotamer_stacks[phi_index, psi_index], rotamer, weight)
            end
        end
    end

    matrices["HIE"] = copy(matrices["HIS"])

    return matrices
end

function create_residue_library_matrices(::Type{T}, filename::String) where {T <: AbstractFloat}

    matrices = Dict{String, BBD_RotamerLibrary}()
    phi_lower_bound = 0.0
    phi_upper_bound = 0.0
    phi_step        = 0.0
    psi_lower_bound = 0.0
    psi_upper_bound = 0.0
    psi_step        = 0.0

    open(filename) do file
        for line in readlines(file)
            if startswith(line, "# phi interval, deg")
                elem = split(line)
                phi_lower_bound = deg2rad(parse(T, elem[5][2:(end-1)]))
                phi_upper_bound = deg2rad(parse(T, elem[6][1:(end-1)]))
            end

            if startswith(line, "# phi step, deg")
                elem = split(line)
                phi_step = deg2rad(parse(T, elem[5]))
            end

            if startswith(line, "# psi interval, deg")
                elem = split(line)
                psi_lower_bound = deg2rad(parse(T, elem[5][2:(end-1)]))
                psi_upper_bound = deg2rad(parse(T, elem[6][1:(end-1)]))
            end

            if startswith(line, "# psi step, deg")
                elem = split(line)
                psi_step = deg2rad(parse(T, elem[5]))
            end

            if startswith(line, "# Input data taken from")
                name = string(split(line)[6])
                phis = collect(phi_lower_bound:phi_step:phi_upper_bound)
                phis[end] = phi_upper_bound
                psis = collect(psi_lower_bound:psi_step:psi_upper_bound)
                psis[end] = psi_upper_bound
                n_phis, n_psis = length(phis), length(psis)
                rotamers = Matrix{RotamerStack}(undef, n_phis, n_psis)
                matrices[name] = BBD_RotamerLibrary(name, phis, psis, rotamers)
            end
        end
    end

    return matrices
end

load_dunbrack(filename::String) = begin
    load_dunbrack(ProtoSyn.Units.defaultFloat, filename)
end
load_dunbrack(::Type{T}) where {T <: AbstractFloat} = begin
    load_dunbrack(T, ProtoSyn.resource_dir * "/Peptides/dunbrack_rotamers.lib")
end
load_dunbrack() = load_dunbrack(ProtoSyn.Units.defaultFloat)