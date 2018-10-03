@doc raw"""
    load_from_gro(i_file::String)::Common.State

Return a new [`Common.State`](@ref) by loading the atom positions and names from the input .gro file.
As a default, `state.energy` is [`NullEnergy`](@ref) and `state.forces` are set to zero.

# Examples
```julia-repl
julia> Common.load_from_gro("molecule.gro")
Common.State(size=2, energy=Null, xyz=[1.1 1.1 1.1; 2.2 2.2 2.2], forces=[0.0 0.0 0.0; 0.0 0.0 0.0], atnames=["C", "O"])
```
See also: [`load_from_pdb`](@ref)
"""
function load_from_gro(i_file::String)::Common.State

    #Initialize empty arrays
    xyz     = Array{Array{Float64, 2}, 1}()
    atnames = Array{String, 1}()

    #Read file (from file name) and recover XYZ and ATOM NAMES
    open(i_file, "r") do f
        for (index, line) in enumerate(eachline(f))
            elem = split(line)
            if length(elem) > 3 && index > 2
                push!(xyz, map(x -> parse(Float64, x), [elem[4] elem[5] elem[6]]))
                push!(atnames, elem[2])
            end
        end
    end

    n = length(xyz)
    return Common.State(n, Common.NullEnergy(), vcat(xyz...), zeros(n, 3), atnames)
end

# ----------------------------------------------------------------------------------------------------------

@doc raw"""
    load_from_pdb(i_file::String)::Common.State

Return a new [`Common.State`](@ref) by loading the atom positions and names from the input .pdb file.
As a default, `state.energy` is [`NullEnergy`](@ref) and `state.forces` are set to zero.

# Examples
```julia-repl
julia> Common.load_from_pdb("molecule.pdb")
Common.State(size=2, energy=Null, xyz=[1.1 1.1 1.1; 2.2 2.2 2.2], forces=[0.0 0.0 0.0; 0.0 0.0 0.0], atnames=["C", "O"])
```
See also: [`load_from_gro`](@ref)
"""
function load_from_pdb(i_file::String)::Common.State

    xyz     = Array{Array{Float64, 2}, 1}()
    atnames = Array{String, 1}()

    open(i_file, "r") do f
        for (index, line) in enumerate(eachline(f))
            elem = split(line)
            if elem[1] == "ATOM"
                push!(xyz, map(x -> parse(Float64, x)/10, [elem[7] elem[8] elem[9]]))
                push!(atnames, elem[3])
            end
        end
    end

    n = length(xyz)
    return Common.State(n, Common.NullEnergy(), vcat(xyz...), zeros(n, 3), atnames)
end