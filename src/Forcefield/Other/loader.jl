function load_contact_maps_from_file(input_file::String, metadata::Common.Metadata)

    #Gather residue information
    resnum2ca::Dict{Int64, Int64} = Dict{Int64, Int64}()
    for residue in metadata.residues
        ca_index = filter(atom -> metadata.atoms[atom].name == "CA", residue.atoms)[1]
        resnum2ca[metadata.atoms[ca_index].res_num] = ca_index
    end

    contact_pairs::Vector{ContactPair} = Vector{ContactPair}()
    open(input_file) do f
        for line in eachline(f)
            elem = split(line)
            if length(elem) > 0 && elem[1] != "i"
                push!(contact_pairs, ContactPair(resnum2ca[parse(Int64, elem[1])], resnum2ca[parse(Int64, elem[2])], parse(Float64, elem[5])))
            end
        end
    end
    return contact_pairs
end