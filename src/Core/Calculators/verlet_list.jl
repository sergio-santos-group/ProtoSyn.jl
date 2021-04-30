"""
    @reduce(n::Int, op::Symbol, ex::Expr)

        Wrapping macro for @ncall. Is equivalent and generates
        op(ex_1, ..., ex_n)

    Example:

        @reduce 2 (+) u -> v_u^2
        > Generates +(v_1^2, v_2^2)
"""
macro reduce(n::Int, op::Symbol, ex::Expr)
    esc(:(@ncall($n, $op, $ex)))
end


"""
    VerletList(size::Int, capacity::Int, cutoff::Float64, offset::Vector{Int}, list::Vector{Int})

Hold information regarding the neighbouring particles of each particle
in the system (with number of particles `size`). A [`VerletList`](@ref) is
comprised of two arrays: the `list`, containing the neighbouring atoms `j` of
atom `i`, with sequential atoms divided by an invalid `-1` entry; and the
`offset` array, containing the positions in array `list` for each atom `i`.
Neighbouring atoms are defined as having a distance bellow the defined `cutoff`.
The main objective of a [`VerletList`](@ref) is to speed up calculations (by
ignoring long-range interactions between [`Atom`](@ref) instances) and to lower
the total amount of memory allocated (the number of allocated [`Atom`](@ref)
entries is at most the `capacity` of the [`VerletList`](@ref)). Note that, given
the motion of particles in a simulation, a [`VerletList`](@ref) can quickly
become obsolete, and needs to be updated using [`update!`](@ref).

    VerletList(size::Int)
    
Creates a new [`VerletList`](@ref) with infinite `cutoff` (holds all atoms in
the molecule).

# Fields
* `size::Int` - The number of [`Atom`](@ref) instances this [`VerletList`](@ref) makes reference to. Should be the size of `:offset` field;
* `capacity::Int` - Maximum number of interaction pairs listed in this [`VerletList`](@ref); 
* `cutoff::Float64` - Interactions are considered when the distance between two [`Atom`](@ref) instances is less than this value;
* `offset::Vector{Int}` - Vector with the starting index for the neighbouring [`Atom`](@ref)`.id` entries in the `:list` field;
* `list::Vector{Int}` - Vector with the neighbouring [`Atom`](@ref)`.id` entries, in sectors separated by invalid entries (such as `-1`).

# See also
[`update!`](@ref) [`distance`](@ref)

# Examples
```jldoctest
julia> ProtoSyn.Calculators.VerletList(4)
ProtoSyn.Calculators.VerletList(4, 4, Inf, [0, 0, 0, 0], [0, 0, 0, 0])
```
"""
mutable struct VerletList
    size::Int              # number of particles the list refers to 
    capacity::Int          # capacity of list 
    cutoff::Float64        # cutoff (in Angstrom Å)
    offset::Vector{Int}    # pointer particle -> list position
    list::Vector{Int}      # nblist
end

VerletList(size::Int) = begin
    return VerletList(size, size, Inf, zeros(Int, size), zeros(Int, size))
end


"""
    resize!(verlet_list::VerletList, n::Int)

        Resize a Verlet list to hold 'n' particles. If 'n' is higher than the
        current length of the list, appends 0's.

    Example:

        resize!(vl, 8)
"""
Base.resize!(verlet_list::VerletList, n::Int) = begin
    resize!(verlet_list.list, n)
    verlet_list.capacity = n
end


"""
    update!([::Type{ProtoSyn.SISD_0}], verlet_list::VerletList, pose::Pose, [selection::Opt{ProtoSyn.AbstractSelection} = nothing])
    update!([::Type{ProtoSyn.SIMD_1}], verlet_list::VerletList, pose::Pose, [selection::Opt{ProtoSyn.AbstractSelection} = nothing])

Updates the given [`VerletList`](@ref) (using a `SISD_0` or `SIMD_1`
acceleration approach) according to the defined 'verlet_list.cutoff' and the
given coordinates in the [`Pose`](@ref) `pose` (in AoS format). If the
acceletarion type is not given, the default `ProtoSyn.acceleration.active` is
employed. If an `AbstractSelection` `selection` is provided, only include the
pair of [`Atom`](@ref) instances as interacting, in the [`VerletList`](@ref), if
both instances are selected. 

# Examples
```jldoctest
julia> update!(verlet_list, pose)
    ...

julia> update!(verlet_list, pose, an"CA")
    ...
```
"""
function update!(::Type{ProtoSyn.SISD_0}, verlet_list::VerletList, pose::Pose, selection::Opt{ProtoSyn.AbstractSelection} = nothing)
    # coords must be in AoS format
    
    @assert verlet_list.size == pose.state.size "incompatible sizes"

    
    # cutoff squared
    cutsq  = convert(eltype(pose.state), verlet_list.cutoff * verlet_list.cutoff)
    offset = 1
    natoms = pose.state.size
    coords = pose.state.x.coords

    if selection !== nothing
        mask = selection(pose)
    else
        mask = !ProtoSyn.Mask{Atom}(natoms)
    end

    @inbounds for i = 1:natoms

        !(mask[i]) && continue # If not selected, continue

        verlet_list.offset[i] = offset
        @nexprs 3 u -> xi_u = coords[u, i]

        for j = (i+1):natoms

            !(mask[j]) && continue # If not selected, continue

            @nexprs 3 u -> vij_u = coords[u, j] - xi_u
            dij_sq = @reduce 3 (+) u -> vij_u * vij_u

            # println("$i - $j: $dij_sq < $cutsq = $(dij_sq < cutsq)")
            if dij_sq < cutsq
                verlet_list.list[offset] = j
                offset += 1

                if offset == verlet_list.capacity
                    resize!(verlet_list, verlet_list.capacity + natoms)
                end
            end
        end

        verlet_list.list[offset] = -1
        offset += 1

        if (i < natoms) && (offset == verlet_list.capacity)
            resize!(verlet_list, verlet_list.capacity + natoms)
        end
    end
    
    return verlet_list
end


function update!(::Type{ProtoSyn.SIMD_1}, verlet_list::VerletList, pose::Pose, selection::Opt{ProtoSyn.AbstractSelection} = nothing)

    @assert verlet_list.size == pose.state.size "incompatible sizes"

    # cutoff squared
    T               = eltype(pose.state)
    cutsq           = convert(T, verlet_list.cutoff * verlet_list.cutoff)
    offset          = 1
    natoms          = pose.state.size
    remaining_mask  = Vec{4, T}((1, 1, 1, 0))

    # For the coords we need to add an extra 0, because the last coord will
    # be loaded in a 4-wide vector for SIMD calculation. The last value is
    # then ignored.
    coords = push!(pose.state.x[:], T(0))

    if selection !== nothing
        mask = selection(pose)
    else
        mask = !ProtoSyn.Mask{Atom}(natoms)
    end

    @inbounds for i = 1:natoms

        !(mask[i]) && continue # If not selected, continue

        verlet_list.offset[i] = offset

        # Atom number to atom position conversion
        _i = i << 2 - (2 + i)
        vi = vload(Vec{4, T}, coords, _i) # Load XYZ (consecutive)

        @inbounds for j = (i+1):natoms

            !(mask[j]) && continue # If not selected, continue
            
            # Atom number to atom position conversion
            _j = j << 2 - (2 + j)

            # Calculate distance, while ignoring the last value in the
            # 4-wide SIMD vector
            rij = (vload(Vec{4, T}, coords, _j) - vi) * remaining_mask   # xi1, yi1, zi1, ?i1
            dij_sq = sum(rij * rij)

            # println("$i - $j: $dij_sq < $cutsq = $(dij_sq < cutsq)")
            if dij_sq < cutsq
                verlet_list.list[offset] = j
                offset += 1

                if offset == verlet_list.capacity
                    resize!(verlet_list, verlet_list.capacity + natoms)
                end
            end
        end # for j

        verlet_list.list[offset] = -1
        offset += 1

        if (i < natoms) && (offset == verlet_list.capacity)
            resize!(verlet_list, verlet_list.capacity + natoms)
        end
    end # for i

    return verlet_list
end

update!(::Type{ProtoSyn.CUDA_2}, verlet_list::VerletList, pose::Pose, selection::Opt{ProtoSyn.AbstractSelection} = nothing) = begin
    @warn "CUDA_2 acceleration not available for `update!` method. Downgrading to SIMD_1 acceleration ..."
    update!(ProtoSyn.SIMD_1, verlet_list, pose, selection)
end

update!(verlet_list::VerletList, pose::Pose, selection::Opt{ProtoSyn.AbstractSelection} = nothing) = begin
    update!(ProtoSyn.acceleration.active, verlet_list, pose, selection)
end


"""
    neighbours(verlet_list::VerletList, atom_index::Int)

Return a list of the neighbouring [`Atom`](@ref) instances of [`Atom`](@ref)
with `:index` `atom_index`, according to the provided [`VerletList`](@ref)
`verlet_list`. 

# Examples
```jldoctest
julia> ProtoSyn.Calculators.neighbours(vl, 1)
22-element Array{Int64,1}:
  2
  3
  4
  5
  6
  ⋮
 23
 24
 25
 26
 30
```
"""
function neighbours(verlet_list::VerletList, atom_index::Int)
    index = 0
    start = verlet_list.offset[atom_index]
    i = verlet_list.list[start]
    n = [i]
    while true
        index += 1
        i = verlet_list.list[start + index]
        i < 0 && break
        push!(n, i)
    end

    return n
end