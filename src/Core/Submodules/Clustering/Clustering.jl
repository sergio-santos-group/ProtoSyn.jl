module Clustering

    using ProtoSyn
    using ProtoSyn.Units
    using Clustering: ClusteringResult
    using LinearAlgebra: tril
    using ProgressMeter

    """
        rmsd_matrix([T], poses::Vector{Pose})
        rmsd_matrix([T], poses::Vector{Pose}, selection::AbstractSelection)
    
    Given a vector of `Pose`s, calculate the pairwise RMSD (see
    `ProtoSyn.rmsd`) values between all members, after alignment (see
    `ProtoSyn.align!` for more information). If an optional type T <:
    AbstractFloat is given, use this type. Otherwise, will use
    ProtoSyn.Units.defaultFloat. If an `AbstractSelection` (`selection`) is
    given, use that sub-set of atoms for the alignment and RMSD calculation.

    # Examples
    ```jldoctest
    julia> a = ProtoSyn.Clustering.rmsd_matrix(poses, an"CA")
    3×3 Array{Float64,2}:
    0.0      7.18633  9.9984
    7.18633  0.0      7.01093
    9.9984   7.01093  0.0
    ```

    # See also
    `ProtoSyn.align!` `ProtoSyn.rmsd` 
    `ProtoSyn.Clustering.complete_diameter_distances` 
    `ProtoSyn.Clustering.average_diameter_distances`
    """
    function rmsd_matrix(::Type{T}, poses::Vector{Pose}, selection::AbstractSelection) where {T <: AbstractFloat}
        N           = length(poses)
        rmsd_matrix = zeros(T, N, N)

        @showprogress "Measuring RMSD matrix - " for i in 1:N
            for j in (i+1):N
                ProtoSyn.align!(poses[i], poses[j], selection)
                rmsd_matrix[i, j] = ProtoSyn.rmsd(poses[i], poses[j], selection)
                rmsd_matrix[j, i] = ProtoSyn.rmsd(poses[i], poses[j], selection)
            end # for
        end # for

        return rmsd_matrix
    end # function

    rmsd_matrix(poses::Vector{Pose}) = begin
        return rmsd_matrix(Units.defaultFloat, poses, TrueSelection{Atom}())
    end

    rmsd_matrix(::Type{T}, poses::Vector{Pose}) where {T <: AbstractFloat} = begin
        return rmsd_matrix(T, poses, TrueSelection{Atom}())
    end

    rmsd_matrix(poses::Vector{Pose}, selection::AbstractSelection) = begin
        return rmsd_matrix(Units.defaultFloat, poses, selection)
    end

    # --------------------------------------------------------------------------
    # Intra and Inter-cluster distances. For more information, see:
    # https://www.geeksforgeeks.org/ml-intercluster-and-intracluster-distance/

    """
        complete_diameter_distances([T], assignments::Vector{Int}, rm::Matrix{T}) where {T <: AbstractFloat}

    Given an `assignment` vector (see `Clustering` package) and a RMSD matrix
    (see `ProtoSyn.Clustering.rmsd_matrix`), return the intra- and inter-cluster
    complete diameter distances. *Note:* The complete diameter distance is the
    distance between the two most remote objects is the same cluster
    (intra-distance) or in two clusters (inter-distance).

    # Examples
    ```jldoctest
    julia> intra, inter = ProtoSyn.Clustering.complete_diameter_distances(clusters.assignments, rmsd_matrix)
    ([8.611851425604167, 13.379156063323348], [0.0 13.509291152844908; 0.0 0.0])
    ```

    # See also
    `ProtoSyn.Clustering.average_diameter_distances`
    """
    function complete_diameter_distances(::Type{T}, assignments::Vector{Int}, rm::Matrix{T}) where {T <: AbstractFloat}
        K     = length(unique(assignments))
        intra = zeros(T, K)
        inter = zeros(T, K, K)

        for k in 1:K
            ca     = assignments.==k
            if count(ca) > 1
                sub_rm = rm[ca, ca]
                intra[k] = maximum(tril(sub_rm, -1))
            end # if
        
            # Inter cluster distance
            for s in (k+1):K
                sca         = assignments.==s
                sub_rm      = rm[ca, sca]
                inter[k, s] = maximum(sub_rm)
            end # for
        end # for

        return intra, inter
    end # function

    complete_diameter_distances(assignments::Vector{Int}, rm::Matrix{T}) where {T <: AbstractFloat} = begin
        complete_diameter_distances(ProtoSyn.Units.defaultFloat, assignments, rm)
    end


    """
        average_diameter_distances([T], assignments::Vector{Int}, rm::Matrix{T}) where {T <: AbstractFloat}

    Given an `assignment` vector (see `Clustering` package) and a RMSD matrix
    (see `ProtoSyn.Clustering.rmsd_matrix`), return the intra- and inter-cluster
    average diameter distances. *Note:* The average diameter distance is the
    average linkage distance between all the objects is the same cluster
    (intra-distance) or all the objects in two clusters (inter-distance).

    # Examples
    ```jldoctest
    julia> intra, inter = ProtoSyn.Clustering.average_diameter_distances(clusters.assignments, rmsd_matrix)
    ([6.19984080693077, 7.566813059033408], [0.0 8.551225161132669; 0.0 0.0])
    ```

    # See also
    `ProtoSyn.Clustering.complete_diameter_distances`
    """
    function average_diameter_distances(::Type{T}, assignments::Vector{Int}, rm::Matrix{T}) where {T <: AbstractFloat}
        K     = length(unique(assignments))
        intra = zeros(T, K)
        inter = zeros(T, K, K)

        for k in 1:K
            ca = assignments.==k
            if count(ca) > 1
                sub_rm          = rm[ca, ca]
                sub_rm_tril     = tril(sub_rm, -1)
                n               = count(x -> x > 0.0, sub_rm_tril)
                intra[k] = sum(sub_rm_tril)/n
            end # if
        
            # Inter cluster distance
            for s in (k+1):K
                sca         = assignments.==s
                sub_rm      = rm[ca, sca]
                inter[k, s] = sum(sub_rm)/length(sub_rm)
            end # for
        end # for

        return intra, inter
    end # function

    average_diameter_distances(assignments::Vector{Int}, rm::Matrix{T}) where {T <: AbstractFloat} = begin
        average_diameter_distances(ProtoSyn.Units.defaultFloat, assignments, rm)
    end

    # --------------------------------------------------------------------------
    # Dunn index and Davies-Bouldin index. For more information, see:
    # https://www.geeksforgeeks.org/dunn-index-and-db-index-cluster-validity-indices-set-1/
    # https://www.hindawi.com/journals/cin/2015/916240/
    # https://jcheminf.biomedcentral.com/articles/10.1186/s13321-017-0208-0

    """
        dunn_index(intra::Vector{T}, inter::Vector{T}) where {T <: AbstractFloat}

    Given a vector with all intra-distances of all clusters considered and an
    inter-distance matrix between all clusters considered, return the Dunn
    Index. *Note:* The higher the Dunn Index, the better is the clustering.

    # Examples
    ```jldoctest
    julia> ProtoSyn.Clustering.dunn_index(intra, inter)
    1.0097267039046134
    ```

    # See also
    `ProtoSyn.Clustering.davies_bouldin_index`
    """
    function dunn_index(intra::Vector{T}, inter::Matrix{T}) where {T <: AbstractFloat}
        return minimum([x for x in inter if x > 0])/maximum(intra)
    end
    
    """
        davies_bouldin_index(intra::Vector{T}, inter::Vector{T}) where {T <: AbstractFloat}

    Given a vector with all intra-distances of all clusters considered and an
    inter-distance matrix between all clusters considered, return the
    Davies-Bouldin Index. *Note:* The lower the Davies-Bouldin Index, the better
    is the clustering.

    # Examples
    ```jldoctest
    julia> ProtoSyn.Clustering.davies_bouldin_index(intra, inter)
    0.8139215907082012
    ```

    # See also
    `ProtoSyn.Clustering.davies_bouldin_index`
    """
    function davies_bouldin_index(intra::Vector{T}, inter::Matrix{T}) where {T <: AbstractFloat}
        dbi = T(0)
        K   = length(intra)
    
        for i in 1:(K-1)
            D = Vector{T}()
            for j in (i+1):K
                i == j && continue
                push!(D, (intra[i] + intra[j])/inter[i, j])
            end # for
    
            dbi += maximum(D)
        end # for
    
        return dbi /= K
    end # function
end # module