using ProtoSyn.Peptides: Dihedral

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
            #  a) Gathe the name of the Rotamer corresponding residue
            name = string(elem[1])

            #  b) Gather chi values/standar deviation
            chis = Dict{DihedralType, Tuple{T, T}}()
            for index in 1:4
                elem[index + 4] == "0" && continue
                value = parse(T, elem[index + 9])
                sd    = parse(T, elem[index + 13])
                chis[getfield(Dihedral, Symbol("chi$index"))] = (value, sd)
            end

            #  c) Gather probability and create the Rotamer
            probability = parse(T, elem[9])
            rotamer = Rotamer(name, probability, chis)

            # Gather the location to append this rotamer to, in the rotamers
            # matrix
            rl = matrices[name]
            phi_index = first(searchsorted(rl.phis, parse(T, elem[2])))
            psi_index = first(searchsorted(rl.psis, parse(T, elem[3])))
        
            # Append to the found location. If necessary, initiate a new vector
            try
                push!(rl.rotamers[phi_index, psi_index], rotamer)
            catch LoadError
                rl.rotamers[phi_index, psi_index] = Vector{Rotamer}()
                push!(rl.rotamers[phi_index, psi_index], rotamer)
            end
        end
    end

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
                phi_lower_bound = parse(T, elem[5][2:(end-1)])
                phi_upper_bound = parse(T, elem[6][1:(end-1)])
            end

            if startswith(line, "# phi step, deg")
                elem = split(line)
                phi_step = parse(T, elem[5])
            end

            if startswith(line, "# psi interval, deg")
                elem = split(line)
                psi_lower_bound = parse(T, elem[5][2:(end-1)])
                psi_upper_bound = parse(T, elem[6][1:(end-1)])
            end

            if startswith(line, "# psi step, deg")
                elem = split(line)
                psi_step = parse(T, elem[5])
            end

            if startswith(line, "# Input data taken from")
                name = string(split(line)[6])
                phis = collect(phi_lower_bound:phi_step:phi_upper_bound)
                psis = collect(psi_lower_bound:psi_step:psi_upper_bound)
                n_phis, n_psis = length(phis), length(psis)
                rotamers = Matrix{Vector{Rotamer}}(undef, n_phis, n_psis)
                matrices[name] = BBD_RotamerLibrary(name, phis, psis, rotamers)
            end
        end
    end

    return matrices
end

load_dunbrack(filename::String) = load_dunbrack(ProtoSyn.Units.defaultFloat, filename)