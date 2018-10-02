#TO DO:
# 1) DOCUMENT THE FUNCTION
function load_from_gro(i_file::String)

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
    return State(n, Common.NullEnergy(), vcat(xyz...), zeros(n, 3), atnames)
end

#TO DO:
# 1) DOCUMENT THE FUNCTION
function load_from_pdb(i_file::String)

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
    return State(n, Common.NullEnergy(), vcat(xyz...), zeros(n, 3), atnames)
end