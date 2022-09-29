using LinearAlgebra

# * Mask or map

# Main function - CUDA_2
function resolve_mask(A::Type{ProtoSyn.CUDA_2},
    pose::Pose,
    e::Union{Array, CuArray},
    f::Union{Array, CuArray},
    update_forces::Bool,
    mask::Union{ProtoSyn.Mask{C}, Matrix{T}}) where {C <: ProtoSyn.AbstractContainer, T <: AbstractFloat}

    @assert size(e) == size(mask) "The used map $(size(mask)) does not match the size of the current selected coords $(size(e)). Check both the applied selection and mask."
    energy = sum(e .* CuArray(mask))
    
    if update_forces
        f = map(*, f, CuArray(repeat(mask, outer = (1, 1, 3))))
        forces = collect(reshape(sum(f, dims = 1), size(f)[1], 3)')
    else
        forces = nothing
    end

    return energy, forces
end

# Main function - SIMD_1 & SISD_0 (mask)
function resolve_mask(A::Union{Type{ProtoSyn.SISD_0}, Type{ProtoSyn.SIMD_1}},
    pose::Pose,
    e::Union{Array, CuArray},
    f::Union{Array, CuArray},
    update_forces::Bool,
    mask::ProtoSyn.Mask{C}) where {C <: ProtoSyn.AbstractContainer}

    @assert size(e) == size(mask) "The used mask $(size(mask)) does not match the size of the current selected coords $(size(e)). Check both the applied selection and mask."
    energy = sum(e .* mask.content)
    
    if update_forces
        f = map(*, f, repeat(mask, outer = (1, 1, 3)))
        forces = collect(reshape(sum(f, dims = 1), size(f)[1], 3)')
    else
        forces = nothing
    end

    return energy, forces
end

# Main function - SIMD_1 & SISD_0 (matrix)
function resolve_mask(A::Union{Type{ProtoSyn.SISD_0}, Type{ProtoSyn.SIMD_1}},
    pose::Pose,
    e::Union{Array, CuArray},
    f::Union{Array, CuArray},
    update_forces::Bool,
    mask::Matrix{T}) where {T <: AbstractFloat}

    @assert size(e) == size(mask) "The used mask $(size(mask)) does not match the size of the current selected coords $(size(e)). Check both the applied selection and mask."
    energy = sum(e .* mask)
    
    if update_forces
        f = map(*, f, repeat(mask, outer = (1, 1, 3)))
        forces = collect(reshape(sum(f, dims = 1), size(f)[1], 3)')
    else
        forces = nothing
    end

    return energy, forces
end

# * Mask function

# Variant for mask::Function; point to main function above
function resolve_mask(A::Type{<: ProtoSyn.AbstractAccelerationType},
    pose::Pose,
    e::Union{Array, CuArray},
    f::Union{Array, CuArray},
    update_forces::Bool,
    mask::Function)

    resolve_mask(A, pose, e, f, update_forces, mask(pose))
end

# * No mask

# Main function
function resolve_mask(A::Type{<: ProtoSyn.AbstractAccelerationType},
    pose::Pose,
    e::Union{Array, CuArray},
    f::Union{Array, CuArray},
    update_forces::Bool,
    mask::Nothing)

    return sum(e), collect(reshape(sum(f, dims = 1), size(f)[1], 3)')
end


# ------------------------------------------------------------------------------
# * Available potential masks

function get_available_masks(m::Module)
    all_functions  = [string(x) for x in names(m, all=true) if x ∉ (:eval, :include, :proxy) && getproperty(m, x) isa Function && !occursin("#", string(x))]
    mask_functions = [x for x in all_functions if occursin("get", x) & occursin("mask", x) & !occursin("masks", x)]
    push!(mask_functions, "load_map")

    return mask_functions
end

get_available_masks() = get_available_masks(ProtoSyn.Calculators)


"""
    show_available_masks([io::IO], [m::Module])

Prints all available masks for potential restraints Module m (defaults to
`ProtoSyn.Calculators`) to the given `IO` `io` (defaults to `stdout`).
Recursivelly searches any inner Module.

# Examples
```
julia> ProtoSyn.Calculators.show_available_masks()
+-------------------------------------------------------+
| Index | Mask function                                 |
+-------------------------------------------------------+
| 1     | get_bonded_mask                               |
| 2     | get_diagonal_mask                             |
+-------------------------------------------------------+
  └── Consider using the `?` menu to learn more about each mask function.
```
"""
function show_available_masks(io::IO, m::Module)
    mask_functions = get_available_masks(m)

    println(io, "+"*repeat("-", 55)*"+")
    @printf(io, "| %-5s | %-45s |\n", "Index", "Mask function")
    println(io, "+"*repeat("-", 55)*"+")
    for (index, mask_function) in enumerate(mask_functions)
        @printf(io, "| %-5d | %-45s |\n", index, mask_function)
    end
    println(io, "+"*repeat("-", 55)*"+")
    println("  └── Consider using the `?` menu to learn more about each mask function.\n")
end

show_available_masks() = show_available_masks(stdout, ProtoSyn.Calculators)


"""
    get_intra_residue_mask(pose::Pose, [selection::Opt{AbstractSelection}])

For all the [`Atom`](@ref) instances in the provided `AbstractSelection`
`selection` (N), return a 2D N x N [`Mask`](@ref) with all the [`Atom`](@ref)
instances of the given [`Pose`](@ref) `pose` not in the same residue selected. 

!!! ukw "Note:"
    This function is rather heavy and has low performance. If no design effort
    is being made (where the sequence changes), the resulting [`Mask`](@ref)
    from this function can and should be re-used (only calculated once). If, for
    a specific application, the `AbstractSelection` `selection` remains constant
    but the [`Mask`](@ref) needs to be re-calculated (for example, because there
    was a design/mutation step, use the _functor_ resulting from
    [`get_intra_residue_mask`](@ref)).

# See also
[`show_available_masks`](@ref)

# Examples
```
julia> ProtoSyn.Calculators.get_intra_residue_mask(pose, !an"^CA\$|^N\$|^C\$|^H\$|^O\$"r)
ProtoSyn.Mask
 ├── Type: Atom
 ├── Size: (1140, 1140)
 ├── Count: 1279946 / 1299600
 └── Content: [0 0 … 1 1; 0 0 … 1 1; … ; 1 1 … 0 0; 1 1 … 0 0]
```
"""
function get_intra_residue_mask(pose::Pose, selection::Opt{AbstractSelection})
    # ! Note: This function is rather heavy. If no design effort is being made
    # ! (where the sequence changes), the resulting map from this function can
    # ! and should be re-used (only calculated once).

    if selection !== nothing
        sele = ProtoSyn.promote(selection, Atom)
    else
        sele = TrueSelection{Atom}()
    end

    s = sele(pose)
    N = count(s)
    mask = ProtoSyn.Mask{Atom}(N, N)
    for r in eachresidue(pose.graph)
        rm1 = ProtoSyn.promote(SerialSelection{Residue}(r.id, :id), Atom)(pose)
        rm1.content = rm1.content[s.content]
        mask |= ProtoSyn.cross2d(rm1)
    end

    return !mask
end

get_intra_residue_mask(pose::Pose) = get_intra_residue_mask(pose, nothing)


"""
    get_diagonal_mask(pose::Pose, [selection::Opt{AbstractSelection}])

For all the atoms in the provided `AbstractSelection` `selection` (N), return a
2D N x N [`Mask`](@ref) with all the [`Atom`](@ref) instances of the given
[`Pose`](@ref) `pose` not in the natural diagonal selected (i.e. ignores same
[`Atom`](@ref) interaction artifacts).

# See also
[`show_available_masks`](@ref)

# Examples
```
julia> ProtoSyn.Calculators.get_diagonal_mask(pose, an"CA")
ProtoSyn.Mask
 ├── Type: Atom
 ├── Size: (73, 73)
 ├── Count: 5256 / 5329
 └── Content: [0 1 … 1 1; 1 0 … 1 1; … ; 1 1 … 0 1; 1 1 … 1 0]
```
"""
function get_diagonal_mask(pose::Pose, selection::Opt{AbstractSelection})

    if selection !== nothing
        sele = ProtoSyn.promote(selection, Atom)
    else
        sele = TrueSelection{Atom}()
    end

    N = count(sele(pose))
    return !ProtoSyn.Mask{Atom}(BitArray(Matrix{Bool}(I, N, N)))
end

get_diagonal_mask(pose::Pose) = get_diagonal_mask(pose, nothing)


"""
    get_bonded_mask(pose::Pose, [selection::Opt{AbstractSelection}])

For all the atoms in the provided `AbstractSelection` `selection` (N), return a
2D N x N [`Mask`](@ref): for each [`Atom`](@ref) instance of the given
[`Pose`](@ref) `pose` mask out all other bonded [`Atom`](@ref) instance.

# See also
[`show_available_masks`](@ref)

# Examples
```
julia> ProtoSyn.Calculators.get_bonded_mask(pose, an"CA")
ProtoSyn.Mask
 ├── Type: Atom
 ├── Size: (73, 73)
 ├── Count: 5256 / 5329
 └── Content: [0 1 … 1 1; 1 0 … 1 1; … ; 1 1 … 0 1; 1 1 … 1 0]
```
"""
function get_bonded_mask(pose::Pose, selection::Opt{AbstractSelection})
    
    if selection !== nothing
        sele = ProtoSyn.promote(selection, Atom)
    else
        sele = TrueSelection{Atom}()
    end

    atoms = sele(pose, gather = true)
    N = length(atoms)
    m = .!Matrix{Bool}(I, N, N)
    for (i, atom_i) in enumerate(atoms)
        for (j, atom_j) in enumerate(atoms)
            i === j && continue

            if atom_j in atom_i.bonds
                m[i, j] = false
            end
        end
    end

    return ProtoSyn.Mask{Atom}(BitArray(m))
end

get_bonded_mask(pose::Pose) = get_bonded_mask(pose, nothing)


"""
    get_upper_triangular_matrix_mask(pose::Pose, [selection::Opt{AbstractSelection}])

For all the atoms in the provided `AbstractSelection` `selection` (N), return a
2D N x N Matrix{T} with the bottom triangular matrix set to 0.0 (including
diagonal) and upper triangular matrix set to 1.0.

# See also
[`show_available_masks`](@ref)

# Examples
```
julia> ProtoSyn.Calculators.get_upper_triangular_matrix_mask(pose, an"CA")
73×73 Matrix{Float64}:
 (...)
```
"""
function get_upper_triangular_matrix_mask(pose::Pose, selection::Opt{AbstractSelection})

    if selection !== nothing
        sele = ProtoSyn.promote(selection, Atom)
    else
        sele = TrueSelection{Atom}()
    end

    T = eltype(pose.state)
    N = count(sele(pose))
    map = Matrix(UpperTriangular(ones(T, (N, N))))
    map[I(N)] .= zero(T)

    return map
end

get_upper_triangular_matrix_mask(pose::Pose) = begin
    get_upper_triangular_matrix_mask(pose, nothing)
end


"""
    get_upper_triangular_matrix_mask(pose::Pose, [selection::Opt{AbstractSelection}])

For all the atoms in the provided `AbstractSelection` `selection` (N), return a
2D N x N Matrix{T} with the bottom triangular matrix set to 0.0 (including
diagonal) and upper triangular matrix set to -1.0.

# See also
[`show_available_masks`](@ref)

# Examples
```
julia> ProtoSyn.Calculators.get_upper_triangular_matrix_inversed_mask(pose, an"CA")
73×73 Matrix{Float64}:
 (...)
```
"""
function get_upper_triangular_matrix_inversed_mask(pose::Pose, selection::Opt{AbstractSelection})

    if selection !== nothing
        sele = ProtoSyn.promote(selection, Atom)
    else
        sele = TrueSelection{Atom}()
    end

    T = eltype(pose.state)
    N = count(sele(pose))
    map = Matrix(UpperTriangular(ones(T, (N, N)) .* T(-1.0)))
    map[I(N)] .= zero(T)

    return map
end

get_upper_triangular_matrix_inversed_mask(pose::Pose) = begin
    get_upper_triangular_matrix_inversed_mask(pose, nothing)
end


"""
    load_PFRMAT_RR_map([::Type{T}], filename::String) where {T <: AbstractFloat}

Load the map in the `filename` file (i.e. Contact Map). The file should be in
PFRMAT RR format (See: [https://predictioncenter.org/casp13/index.cgi?page=format#RR](https://predictioncenter.org/casp13/index.cgi?page=format#RR)).
Returns an N x N map of the found weights, with pairs not identified in the file
set to 0.0 (N is the maximum indentifier found on the file. As an example, it
might be the case where a peptide has 74 residues, but no pair with residue 74
is found on the file, the maximum identifier found might be 72, for example. In
this case, the resulting map will have size 72 x 72. In order to ensure the
loaded map size matches the underlying peptide size, consider adding an entry of
0.0 on the map file, with the correct maximum identifier). Note: If no
optional type `T` is provided, will use `ProtoSyn.Units.defaultFloat`.

# Examples
```
julia> cmap = ProtoSyn.Calculators.load_PFRMAT_RR_map("contact_map_example.txt")
73×73 Array{Float64,2}:
 (...)
```
"""
function load_PFRMAT_RR_map(::Type{T}, filename::String, i::Int) where {T <: AbstractFloat}

    _map = Dict{Tuple{Int, Int}, T}()

    open(filename, "r") do map_file
        for line in eachline(map_file)
            elems = split(line)
            tryparse(Int, elems[1]) === nothing && continue

            α = parse(T, elems[i])
            _map[(parse(Int, elems[1]), parse(Int, elems[2]))] = α
        end
    end

    k    = keys(_map)
    N    = maximum((maximum(last, k), maximum(first, k)))
    map  = zeros((N, N))

    for (key, value) in _map
        map[key[1], key[2]] = value
        map[key[2], key[1]] = value
    end

    return map
end

load_PFRMAT_RR_map(filename::String, i::Int) = load_PFRMAT_RR_map(ProtoSyn.Units.defaultFloat, filename, i)
load_PFRMAT_RR_map(filename::String)         = load_PFRMAT_RR_map(ProtoSyn.Units.defaultFloat, filename, 3)