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

"""
    intra_residue_mask(pose::Pose, selection::AbstractSelection)

For all the atoms in the provided `AbstractSelection` `selection` (N), return a
2D N x N [`Mask`](@ref) with all the atoms of the given [`Pose`](@ref) `pose`
not in the same residue selected. 

!!! ukw "Note:"
    This function is rather heavy and has low performance. If no design effort
    is being made (where the sequence changes), the resulting [`Mask`](@ref)
    from this function can and should be re-used (only calculated once). If, for
    a specific application, the `AbstractSelection` `selection` remains constant
    but the [`Mask`](@ref) needs to be re-calculated (for example, because there
    was a design/mutation step, use the _functor_ resulting from
    [`get_intra_residue_mask`](@ref)).

# See also
[`diagonal_mask`](@ref) [`get_intra_residue_mask`](@ref)

# Examples
```
julia> ProtoSyn.Calculators.intra_residue_mask(pose, !an"^CA\$|^N\$|^C\$|^H\$|^O\$"r)
ProtoSyn.Mask{Atom}(Bool[1 1 … 0 0; 1 1 … 0 0; … ; 0 0 … 1 1; 0 0 … 1 1])
```
"""
function intra_residue_mask(pose::Pose, selection::AbstractSelection)
    # ! Note: This function is rather heavy. If no design effort is being made
    # ! (where the sequence changes), the resulting map from this function can
    # ! and should be re-used (only calculated once).

    s = selection(pose)
    N = count(s)
    mask = ProtoSyn.Mask{Atom}(N, N)
    for r in eachresidue(pose.graph)
        rm1 = ProtoSyn.promote(SerialSelection{Residue}(r.id, :id), Atom)(pose)
        rm1.content = rm1.content[s.content]
        mask |= ProtoSyn.cross2d(rm1)
    end

    return !mask
end


"""
    get_intra_residue_mask(selection::AbstractSelection)

Provides the [`intra_residue_mask`](@ref) function as a _functor_, which will
calculate the intra residue mask for the given `AbstractSelection` `selection`.
(Only the [`Atom`](@ref) instances in the selection are considered, all other
atoms are not included in the [`Mask`](@ref)). Useful when creating a new
[`EnergyFunctionComponent`](@ref) or when the [`Mask`](@ref) should be updated
each step/call.

# See also
[`intra_residue_mask`](@ref)

# Examples
```jldoctest
julia> ProtoSyn.Calculators.get_intra_residue_mask(!an"^CA\$|^N\$|^C\$|^H\$|^O\$"r)
(::ProtoSyn.Calculators.var"#_intra_residue_mask#5"{UnarySelection{ProtoSyn.Stateless}}) (generic function with 1 method)
```
"""
function get_intra_residue_mask(selection::AbstractSelection)
    return function _intra_residue_mask(pose::Pose)
        s = selection(pose)
        N = count(s)
        mask = ProtoSyn.Mask{Atom}(N, N)
        for r in eachresidue(pose.graph)
            rm1 = ProtoSyn.promote(SerialSelection{Residue}(r.id, :id), Atom)(pose)
            rm1.content = rm1.content[s.content]
            mask |= ProtoSyn.cross2d(rm1)
        end

        return !mask
    end
end


"""
    diagonal_mask(pose::Pose, selection::AbstractSelection)

For all the atoms in the provided `AbstractSelection` `selection` (N), return a
2D N x N [`Mask`](@ref) with all the [`Atom`](@ref) instances of the given
[`Pose`](@ref) `pose` not in the natural diagonal selected (i.e. ignores same
atom interaction artifacts).

!!! ukw "Note:"
    When the selection is constant but the resulting [`Mask`](@ref) needs to be
    re-calculated every call/step (for example, due to a design or mutation
    step), consider using the _functor_ from [`get_diagonal_mask`](@ref).

# See also
[`intra_residue_mask`](@ref) [`get_diagonal_mask`](@ref)

# Examples
```
julia> ProtoSyn.Calculators.diagonal_mask(pose, an"CA")
ProtoSyn.Mask{Atom}(3, 3)
3×3 BitArray{2}:
 0  1  1
 1  0  1
 1  1  0
```
"""
function diagonal_mask(pose::Pose, selection::AbstractSelection)

    N = count(selection(pose))
    return !ProtoSyn.Mask{Atom}(BitArray(Matrix{Bool}(I, N, N)))
end


"""
    get_diagonal_mask(selection::AbstractSelection)

Provides the [`diagonal_mask`](@ref) as a _functor_, which will calculate the
diagonal mask for the given `AbstractSelection` `selection`. Useful when
creating a new [`EnergyFunctionComponent`](@ref) or when the [`Mask`](@ref)
should be updated each step/call.

# See also
[`diagonal_mask`](@ref)

# Examples
```
julia> ProtoSyn.Calculators.get_diagonal_mask(an"CA")
(::ProtoSyn.Calculators.var"#_diagonal_mask#6"{FieldSelection{ProtoSyn.Stateless,Atom}}) (generic function with 1 method)
```
"""
function get_diagonal_mask(selection::AbstractSelection)
    return function _diagonal_mask(pose::Pose)
        N = count(selection(pose))
        return !ProtoSyn.Mask{Atom}(BitArray(Matrix{Bool}(I, N, N)))
    end
end

function get_diagonal_mask()
    return function _diagonal_mask(pose::Pose)
        N = pose.state.size
        return !ProtoSyn.Mask{Atom}(BitArray(Matrix{Bool}(I, N, N)))
    end
end

# ---
# DEV
function get_bonded_mask()
    return function _diagonal_mask(pose::Pose)
        N = pose.state.size
        m = .!Matrix{Bool}(I, N, N)
        for atom in eachatom(pose.graph)
            i = atom.index
            for bond in atom.bonds
                j = bond.index
                m[i, j] = false
            end
        end

        return ProtoSyn.Mask{Atom}(BitArray(m))
    end
end


"""
    load_map([::Type{T}], filename::String) where {T <: AbstractFloat}

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
julia> cmap = ProtoSyn.Calculators.load_map("contact_map_example.txt")
73×73 Array{Float64,2}:
 ...
```
"""
function load_map(::Type{T}, filename::String) where {T <: AbstractFloat}

    _map = Dict{Tuple{Int, Int}, T}()

    open(filename, "r") do map_file
        for line in eachline(map_file)
            elems = split(line)
            length(elems[1]) > 4 && continue
            elems[1] == "END" && continue
            try
                α = parse(T, elems[5])
                _map[(parse(Int, elems[1]), parse(Int, elems[2]))] = α
            catch BoundsError
                continue
            end
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

load_map(filename::String) = load_map(ProtoSyn.Units.defaultFloat, filename)