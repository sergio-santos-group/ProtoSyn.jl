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
    printstyled("(SETUP) ▲ Loaded $(length(solv_pairs)) solvation pairs\n", color = 9)
    return solv_pairs
end

@doc raw"""
    compile_solv_pairs(phi_dihedrals::Vector{Common.Dihedral}, λ_eSol::Float64)::Vector{SolvPair}

Apply default solvation coeficients to the correspondent residue in [`SolvPair`] (@ref)'s.

# Examples
```julia-repl
julia> Forcefield.CoarseGrain.compile_solv_pairs(phi_dihedrals, 0.01)
[Forcefield.CoarseGrain.SolvPair(i=1, coef=-3.5), Forcefield.CoarseGrain.SolvPair(i=2, coef=2.0), ...]
```
"""
function compile_solv_pairs(phi_dihedrals::Vector{Common.Dihedral}, λ_eSol::Float64)::Vector{SolvPair}
    solv_pairs = map(x -> SolvPair(x.a3, default_aa_coef[string(x.residue.name[1])] * λ_eSol), phi_dihedrals)
    printstyled("(SETUP) ▲ Compiled $(length(solv_pairs)) solvation pairs\n", color = 9)
    return solv_pairs
end


@doc raw"""
    compile_hb_groups(atoms::Vector{Common.AtomMetadata}, λ_eSol::Float64)::Vector{HbGroup}

Compile hydrogen bonding groups from atom metadata.

# Examples
```julia-repl
julia> Forcefield.CoarseGrain.compile_hb_groups(metadata.atoms, 1.0)
[Forcefield.CoarseGrain.HbGroup(n=1, h=2, c=4, o=5, coef=1.0), Forcefield.CoarseGrain.HbGroup(n=6, h=7, c=9, o=10, coef=1.0), ...]
```
"""
function compile_hb_groups(atoms::Vector{Common.AtomMetadata}, λ_eSol::Float64)::Vector{HbGroup}

    hb_groups::Vector{HbGroup} = HbGroup[]
    for residue in Common.iterate_by_residue(atoms)
        if Aux.conv123(residue[1].res_name) == "PRO"
            continue
        end
        push!(hb_groups, HbGroup(
            filter(atom -> atom.name == "N", residue)[1].index,
            filter(atom -> atom.name == "H", residue)[1].index,
            filter(atom -> atom.name == "C", residue)[1].index,
            filter(atom -> atom.name == "O", residue)[1].index,
            λ_eSol))
    end
    printstyled("(SETUP) ▲ Compiled $(length(hb_groups)) hydrogen bonding groups\n", color = 9)
    return hb_groups
end