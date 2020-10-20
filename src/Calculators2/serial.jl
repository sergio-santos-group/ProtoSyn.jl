function serial(coords::Matrix{T}) where {T <: AbstractFloat}
    # coords must be in AoS format

    natoms = size(coords)[2]
    distance_matrix = zeros(natoms, natoms)

    for i in 1:(natoms - 1)

        for j in (i + 1):natoms
            rij = coords[:, j] - coords[:, i]
            d = sqrt(sum(rij.*rij))
            distance_matrix[i, j] = d
        end
    end
    
    return distance_matrix
end

serial(pose::Pose) = serial(pose.state.x)
serial(state::State{T}) where {T <: AbstractFloat} = serial(state.x)

serial(pose::Pose, selection::ProtoSyn.AbstractSelection) = begin

    s = selection(pose, gather = true)
    coords = zeros(eltype(pose.state), 3, length(s))
    for (i, atom) in enumerate(s)
        coords[:, i] = pose.state[atom].t
    end

    return serial(coords)
end


function serial(state::State{T}, vlist::VerletList) where {T <: AbstractFloat}
    
    natoms = state.size
    coords = state.x
    distance_matrix = zeros(natoms, natoms)
    
    @inbounds for i = 1:natoms-1

        # ptr -> location of the first neighbor of atom i
        # ptr_stop -> location of the last neighbor of atom i
        ptr = vlist.offset[i]
        if vlist.list[ptr] < 1
            continue
        end
        ptr_stop = vlist.offset[i + 1] - 2

        # load coordinates for the i-th atom
        @nexprs 3 u -> ri_u = coords[u, i]
        
        while ptr <= ptr_stop
            j = vlist.list[ptr]

            # load coordinates for the j-th atom
            # and calculate the ij vector
            @nexprs 3 u -> rij_u = coords[u, j] - ri_u
            
            # calculate the squared distance. Skip
            # if greater than cutoff
            dij_sq = @reduce 3 (+) u -> rij_u*rij_u
            distance_matrix[i, j] = sqrt(dij_sq)
            
            ptr += 1
        end
    end
    
    return distance_matrix
end