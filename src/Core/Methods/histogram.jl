"""
# TODO DOCUMENTATION
# USAGE IN GENERALIZED BORN NN
"""
function hist_by_distance_by_elem(pose::Pose, selection::Opt{AbstractSelection} = nothing; cutoff::T = 16.0, bin::T = 0.2, elems::Vector{String} = ["H", "O", "N", "C", "S"], dm::Opt{Matrix{T}} = nothing) where {T <: AbstractFloat}

    if selection === nothing
        sele = TrueSelection{Atom}()
    else
        sele = ProtoSyn.promote(selection, Atom)
    end

    # Pre-calculate distance matrix
    if dm === nothing
        dm = collect(ProtoSyn.Calculators.distance_matrix(pose, selection))
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
            init_elem_pos = elempos[atoms[j].symbol]

            # Find the bin position after the initial element position
            bin_pos = ceil(Int, d / bin)

            # Increment count of this element in this distance bin
            hist[i, init_elem_pos + bin_pos] += 1
        end
    end

    return hist
end