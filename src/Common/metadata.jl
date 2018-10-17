# ----------------------------------------------------------------------------------------------------------
#                                                METADATA

@doc raw"""
    AtomMetadata(name::String[, elem::String = name, res_num::Int64 = 1, res_name::String = "UNK", chain_id::Union{String, Nothing} = nothing, connects::Union{Vector{Int64}, Nothing} = nothing])

Define an atom metadata, containing extra information pertaining the [`State`](@ref).

# Arguments
- `name::String`: Name of the atom.
- `elem::String`: (Optional) Element of the atom (Default: `name`).
- `res_num::Int64`: (Optional) Number of the residue this atom belongs to (Default: 1).
- `res_name::Union{String, Nothing}`: (Optional) Name of the residue this atom belongs to (Default: "UNK").
- `chain_id::String`: (Optional) Name of the chain that contains the residue this atom belongs to (Default: nothing).
- `connects::Union{Vector{Int64}, Nothing}`: (Optional) List of *global* atom indices that this atom is connected to (Default: nothing). 

# Examples
```julia-repl
julia> AtomMetadata("H1", "H", 2, "VAL", "A", [4])
AtomMetadata(name=H1, elem=H, res_num=2, res_name=VAL, chain_id=A, connects=[4])

julia> AtomMetadata("H1")
AtomMetadata(name=H1, elem=H1, res_num=1, res_name=UNK, chain_id=nothing, connects=nothing)
```
See also: [`iter`](@ref)
"""
mutable struct AtomMetadata

    name::String
    elem::String
    res_num::Int64
    res_name::String
    chain_id::Union{String, Nothing}
    connects::Union{Vector{Int64}, Nothing}

    AtomMetadata(name::String; elem::String = name, res_num::Int64 = 1, res_name::String = "UNK", chain_id::Union{String, Nothing} = nothing, connects::Union{Vector{Int64}, Nothing} = nothing) = new(name, elem, res_num, res_name, chain_id, connects)
end
function Base.show(io::IO, b::AtomMetadata) #SHOULD BE IMPROVED
    if b.chain_id != nothing && b.connects != nothing
        print(io, "AtomMetadata(name=$(b.name), elem=$(b.elem), res_num=$(b.res_num), res_name=$(b.res_name), chain_id=$(b.chain_id), connects=$(b.connects))")
    elseif b.chain_id != nothing
        print(io, "AtomMetadata(name=$(b.name), elem=$(b.elem), res_num=$(b.res_num), res_name=$(b.res_name), chain_id=$(b.chain_id), connects=nothing)")
    elseif b.connects != nothing
        print(io, "AtomMetadata(name=$(b.name), elem=$(b.elem), res_num=$(b.res_num), res_name=$(b.res_name), chain_id=nothing, connects=$(b.connects))")
    else
        print(io, "AtomMetadata(name=$(b.name), elem=$(b.elem), res_num=$(b.res_num), res_name=$(b.res_name), chain_id=nothing, connects=nothing)")
    end
end
function Array{AtomMetadata, 1}(n::Int64)
    o = Vector{AtomMetadata}()
    for i in 1:n
        push!(o, AtomMetadata("_"))
    end
    return o
end


@doc raw"""
    iter(data::Vector{AtomMetadata}; property::Symbol = :res_num)

Iterate over an array of AtomMetadata objects, grouping them based on `property` (Default: `:res_num`)

# Examples
```julia-repl
julia> for residue in iter(state.metadata)
    println(residue)
end
[AtomMetadata(name=H1, elem=H1, res_num=1, res_name=UNK, chain_id=nothing, connects=nothing), AtomMetadata(name=H2, elem=H2, res_num=1, res_name=UNK, chain_id=nothing, connects=nothing)]
[AtomMetadata(name=H3, elem=H3, res_num=2, res_name=UNK, chain_id=nothing, connects=nothing), AtomMetadata(name=H4, elem=H4, res_num=2, res_name=UNK, chain_id=nothing, connects=nothing)]
```
See also: [`AtomMetadata`](@ref)
"""
function iter(data::Vector{AtomMetadata}; property::Symbol = :res_num)

    #GET UNIQUE VALUES FOR THIS PROPERTY
    up = Vector{Any}()
    for s in data
        getproperty(s, property) in up ? nothing : push!(up, getproperty(s, property))
    end

    #GROUP ELEMENTS BASED ON THE UNIQUE VALUES FOR THIS PROPERTY
    f = Vector{Vector{AtomMetadata}}()
    idx::Int64 = 1
    cmp::Any = up[idx]
    push!(f, Vector{AtomMetadata}())
    for s in data
        if getproperty(s, property) == cmp
            push!(f[idx], s)
        else
            push!(f, Vector{AtomMetadata}())
            idx += 1
            cmp = up[idx]
            push!(f[idx], s)
        end
    end
    return f
end