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
in the system. A Verlet list is comprised of two arrays: the `list`, containing
the neighbouring atoms (j) of atom i, with sequential atoms divided by a `-1`
entry; and the `offset` array, containing the positions in array `list` for each
atom i. Neighbouring atoms are defined as having a distance bellow the defined
`cutoff`.

    VerletList(size::Int)
    
Creates a new Verlet list with infinite cutoff (holds all atoms in the
molecule).

# See also
`update_serial!`, `update_simd!`

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

VerletList(size::Int) = VerletList(size, size, Inf, zeros(Int, size), zeros(Int, size))


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
    verlet_list
end


"""
    update_serial!(verlet_list::VerletList, pose::Pose, selection::Opt{ProtoSyn.AbstractSelection} = nothing)

Updates the Verlet list (using a serial SISD approach) according to the
defined 'verlet_list.cutoff' and the given coordinates 'coords' (must be
in AoS format).

# Examples
```jldoctest
julia> update_serial!(verlet_list, pose.state.x)
```
"""
function update_serial!(verlet_list::VerletList, pose::Pose, selection::Opt{ProtoSyn.AbstractSelection} = nothing)
    # coords must be in AoS format
    
    @assert verlet_list.size == pose.state.size "incompatible sizes"

    
    # cutoff squared
    cutsq  = convert(eltype(pose.state), verlet_list.cutoff * verlet_list.cutoff)
    offset = 1
    natoms = pose.state.size
    coords = pose.state.x

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


export update_simd! 

"""
    update_simd!(verlet_list::VerletList, pose::Pose, selection::Opt{ProtoSyn.AbstractSelection} = nothing) where {T <: AbstractFloat}

Updates the Verlet list (using a serial SISD approach) according to the defined
'verlet_list.cutoff' and the given coordinates defined in `pose.state.x` (must
be in AoS format). If a `selection` is given, only atoms included in that
selection are considered when updating the `verlet_list`.

# Examples
```jldoctest
julia> update_serial!(verlet_list, pose.state.x)
```
"""
function update_simd!(verlet_list::VerletList, pose::Pose, selection::Opt{ProtoSyn.AbstractSelection} = nothing)

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