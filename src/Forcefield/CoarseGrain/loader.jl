@doc raw"""
    load_solv_pairs_from_file(c_αs::Vector{Int64}, λ_eSol::Float64, input_file::String)::Vector{SolvPair}

Read an RSA map file and apply each solvation coeficient to the correspondent residue in [`SolvPair`] (@ref)'s.

# Examples
```julia-repl
julia> Forcefield.CoarseGrain.load_solv_pairs_from_file(c_αs, 0.01, "rsa_map.txt")
[Forcefield.CoarseGrain.SolvPair(i=1, coef=-3.5), Forcefield.CoarseGrain.SolvPair(i=2, coef=2.0), ...]
```
"""
function load_solv_pairs_from_file(c_αs::Vector{Int64}, λ_eSol::Float64, input_file::String)::Vector{SolvPair}

    coefs::Vector{Int64} = Int64[]
    confidences::Vector{Int64} = Int64[]
    solv_pairs::Vector{SolvPair} = SolvPair[]
    open(input_file) do f
        index::Int64 = 1
        for line in eachline(f)
            if index == 1
                coefs = [parse(Int64, x) for x in line]
                coefs .-= 3
            else
                confidences = [parse(Int64, x) for x in line]
            end
            index += 1
        end
    end

    for (index, c_α) in enumerate(c_αs)
        push!(solv_pairs, SolvPair(c_α, coefs[index] * confidences[index] * λ_eSol))
    end
    return solv_pairs
end

@doc raw"""
    load_solv_pairs_default(phi_dihedrals::Vector{Common.Dihedral}, λ_eSol::Float64)::Vector{SolvPair}

Apply default solvation coeficients to the correspondent residue in [`SolvPair`] (@ref)'s.

# Examples
```julia-repl
julia> Forcefield.CoarseGrain.load_solv_pairs_default(phi_dihedrals, 0.01)
[Forcefield.CoarseGrain.SolvPair(i=1, coef=-3.5), Forcefield.CoarseGrain.SolvPair(i=2, coef=2.0), ...]
```
"""
function load_solv_pairs_default(phi_dihedrals::Vector{Common.Dihedral}, λ_eSol::Float64)::Vector{SolvPair}
    return map(x -> SolvPair(x.a3, default_aa_coef[string(x.residue.name[1])] * λ_eSol), phi_dihedrals)
end