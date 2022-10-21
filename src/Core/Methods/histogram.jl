"""
    hist_by_distance_by_elem(pose::Pose, [selection::Opt{AbstractSelection} = nothing]; [cutoff::T = 16.0], [bin::T = 0.2], [elems::Vector{String} = ["H", "O", "N", "C", "S"]], [dm::Opt{Matrix{T}} = nothing]) where {T <: AbstractFloat}

Calculates the given [`Pose`](@ref) `pose` distance histogram by [`Atom`](@ref)
element. The output format matches the requirements for
[`predict_igbr_nn_born_radii`](@ref ProtoSyn.Calculators.GB.predict_igbr_nn_born_radii)
(for more information, see https://academic.oup.com/bioinformatics/article/36/6/1757/5613804).
If an `AbstractSelection` `selection` is provided, only the selected subset of
[`Atom`](@ref) instances is considered for the histogram calculation (promoted
to [`Atom`](@ref) level using the [`promote`](@ref ProtoSyn.promote) method).
The distance histogram is generated for all [`Atom`](@ref) instances within
`cutoff` Å of the focused [`Atom`](@ref), sub-divided in bins of `bin` Å. If a
[`full_distance_matrix`](@ref ProtoSyn.Calculators.full_distance_matrix) `dm` is provided, use it for distance
calculation, otherwise a new distance matrix is generated. Each [`Atom`](@ref)
histogram is divided by elements: one histogram is generated for each atomic
element in `elems` containing all neighbouring [`Atom`](@ref) instances of that
type. 

# Examples
```jldoctest
julia> ProtoSyn.hist_by_distance_by_elem(pose)
1140×400 Matrix{Int64}:
 (...)
```

"""
function hist_by_distance_by_elem(pose::Pose, selection::Opt{AbstractSelection} = nothing; cutoff::T = 16.0, bin::T = 0.2, elems::Vector{String} = ["H", "O", "N", "C", "S"], dm::Opt{Matrix{T}} = nothing) where {T <: AbstractFloat}

    if selection === nothing
        sele = TrueSelection{Atom}()
    else
        sele = ProtoSyn.promote(selection, Atom)
    end

    # Pre-calculate distance matrix
    if dm === nothing
        dm = collect(ProtoSyn.Calculators.full_distance_matrix(pose, sele))
    end

    # Define useful counts
    atoms  = sele(pose, gather = true)
    natoms = length(atoms)
    nbins  = convert(Int, cutoff / bin)
    nelems = length(elems)

    # Pre-set table lookups for initial element position in the histogram
    elempos = Dict{String, Int64}()
    for (i, elem) in enumerate(elems)
        elempos[elem] = nbins * (i - 1)
    end

    # Allocate histogram memory
    hist = zeros(Int64, natoms, nelems * nbins)

    # Loop over each atom pair
    for i in 1:natoms
        for j in 1:natoms

            # Ignore same atoms
            i === j && continue

            # Lookup distance from pre-calculated distance matrix
            d = dm[i, j]

            # Apply cutoff
            d > cutoff && continue

            # Find the initial element position in the histogram
            s = atoms[j].symbol
            if !(s in elems)
                continue
            end
            init_elem_pos = elempos[s]
            
            # Find the bin position after the initial element position
            bin_pos = ceil(Int, d / bin)
            
            if bin_pos === 0 # Case distance between particles is 0.0
                continue
            end

            # Increment count of this element in this distance bin
            hist[i, init_elem_pos + bin_pos] += 1
        end
    end

    return hist
end