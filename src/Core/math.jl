using Base.Cartesian

export rotmat, rotmat!, rotate!, measure

@inline function rotmat!(R::Matrix{Float64}, axis::Vector{Float64}, θ::Float64)
    # q0 = cos(0.5 * θ)
    # q1, q2, q3 = sin(0.5 * θ) * axis ./ norm(axis)
    λ  = sin(0.5 * θ) / norm(axis)
    q0 = cos(0.5 * θ)
    q1 = λ * axis[1]
    q2 = λ * axis[2]
    q3 = λ * axis[3]
    R[1,1] = 1-2*q2*q2-2*q3*q3
    R[1,2] =   2*q1*q2-2*q0*q3
    R[1,3] =   2*q1*q3+2*q0*q2
    R[2,1] =   2*q2*q1+2*q0*q3
    R[2,2] = 1-2*q3*q3-2*q1*q1
    R[2,3] =   2*q2*q3-2*q0*q1
    R[3,1] =   2*q3*q1-2*q0*q2
    R[3,2] =   2*q3*q2+2*q0*q1
    R[3,3] = 1-2*q1*q1-2*q2*q2
    R
end

@inline rotmat(axis::Vector{Float64}, θ::Float64) = rotmat!(zeros(3,3), axis, θ)



# function kabsch(reference::Matrix{Float64}, mobile::Matrix{Float64})
#     # shift to origin both matrices
#     refmean = mean(reference, dims=1)
#     mobmean = mean(mobile, dims=1)
#     reference .-= refmean
#     mobile .-= mobmean
    
#     # calculate correlation matrix
#     corr = mobile' * reference

#     reference .+= refmean
#     # decompose correlation matrix using SVD
#     f = svd(corr)

#     # sign(det(f.V) * det(f.U))
#     if det(f.Vt) * det(f.U) < 0.0
#         f.Vt[3,:] .*= -1.0
#     end

#     # Xr = (R * X')'     <=>
#     # Xr = X * R'        <=>
#     # Xr = X * (V * U')' <=>
#     # Xr = X * (U * V')  <=>
#     # Xr = X *  U * V'

#     return mobile * f.U * f.Vt
# end



# function rotate!(state::State, mol::Molecule,
#     block::AbstractAxisRotatableBlock, θ::Float64)

#     graph = mol.bonds
#     xyz = state.coords

#     axis = [0.0, 0.0, 0.0]

#     # origin of rotation
#     @nexprs 3 k -> o_k = xyz[block.a1,k]
    
#     # axis of rotation along vector a1->a2
#     @nexprs 3 k -> axis[k] = xyz[block.a2,k] - o_k
    
#     # get rotation matrix
#     R = rotmat(axis, θ)
    
#     # perform rotation for all atoms in graph
#     # reachable from "start" using recursion (needs timing tests!)
#     visited = falses(size(xyz,1))
#     function traverse(j::Int)
#         @inbounds begin
#             visited[j] = true
#             @nexprs 3 k -> x_k = xyz[j,k]-o_k
#             @nexprs 3 k -> xyz[j,k] = R[k,1]*x_1 + R[k,2]*x_2 + R[k,3]*x_3 + o_k
#             for i in graph[j]
#                 !visited[i] && traverse(i)
#             end
#         end
#     end
    
#     # make sure to mark a1 & a2 as visited to prevent traversal
#     # to the wrong side before doing the job. It might happen that
#     # start is equal to a2 in which case, despited being already marked
#     # as visited, traversal will still take place because rotation 
#     # is performed on the current index (j in function traverse) without
#     # any check.
#     visited[block.a1] = true
#     visited[block.a2] = true
#     traverse(block.start)

#     # return coordinates
#     state
# end

rotate!(state::State, mol::Molecule, block::AbstractAxisRotatableBlock, θ::Float64) = begin
    rotate!(state.coords, mol.bonds, block, θ, falses(mol.size), zeros(3,3))
    state
end

# rotate!(xyz::Matrix{Float64}, mol::Molecule, block::AbstractAxisRotatableBlock, θ::Float64) = begin
#     rotate!(xyz, mol, block, θ, falses(mol.size), zeros(3,3))
# end

function rotate!(xyz::Matrix{Float64}, graph::ConnectGraph,
    block::AbstractAxisRotatableBlock, θ::Float64, visited::BitVector, rmat::Matrix{Float64})

    @assert length(visited) == length(graph) "bitvec and mol must have same size"

    axis = [0.0, 0.0, 0.0]

    # origin of rotation
    @nexprs 3 k -> o_k = xyz[block.a1, k]
    
    # axis of rotation along vector a1->a2
    @nexprs 3 k -> axis[k] = xyz[block.a2,k] - o_k
    
    # get rotation matrix
    rotmat!(rmat, axis, θ)
    
    # perform rotation for all atoms in graph
    # reachable from "start" using recursion (needs timing tests!)
    # make sure to mark a1 & a2 as visited to prevent traversal
    # to the wrong side before doing the job. It might happen that
    # start is equal to a2 in which case, despited being already marked
    # as visited, traversal will still take place because rotation 
    # is performed on the current index (j in function traverse) without
    # any check.
    
    fill!(visited, false)
    visited[block.a1] = true
    visited[block.a2] = true
    _rotatetree!(xyz, visited, block.start, rmat, o_1, o_2, o_3, graph)
    
    # return changed state
    xyz
end


function _rotatetree!(xyz::Matrix{Float64}, visited::BitVector, pivot::Int,
    rmat::Matrix{Float64}, o_1::Float64, o_2::Float64, o_3::Float64, graph::ConnectGraph)
    @inbounds begin
        visited[pivot] = true
        @nexprs 3 k -> x_k = xyz[pivot,k] - o_k
        @nexprs 3 k -> xyz[pivot,k] = rmat[k,1]*x_1 + rmat[k,2]*x_2 + rmat[k,3]*x_3 + o_k
    end
    for i in graph[pivot]
        @inbounds !visited[i] && _rotatetree!(xyz, visited, i, rmat, o_1, o_2, o_3, graph)
    end
end


# function rotate!(state::State, mol::Molecule, dihd::Dihedral, θ::Float64)
#     graph = mol.bonds
#     xyz = state.coords

#     axis = [0.0, 0.0, 0.0]

#     # origin of rotation
#     @nexprs 3 k -> o_k = xyz[dihd.a3,k]
    
#     # axis of rotation along vector a2->a3
#     @nexprs 3 k -> axis[k] = o_k - xyz[dihd.a3,k]
    
#     # get rotation matrix
#     R = rotmat(axis, θ)
    
#     # perform rotation for all atoms in graph
#     # reachable from a3 using recursion (needs timing tests!)
#     visited = falses(size(xyz,1))
#     function traverse(j::Int)
#         @inbounds begin
#             visited[j] = true
#             @nexprs 3 k -> x_k = xyz[j,k]-o_k
#             @nexprs 3 k -> xyz[j,k] = R[k,1]*x_1 + R[k,2]*x_2 + R[k,3]*x_3 + o_k
#             for i in graph[j]
#                 !visited[i] && traverse(i)
#             end
#         end
#     end
    
#     # make sure to mark a2 as visited to prevent traversal
#     # to the wrong side before doing the job
#     visited[dihd.a2] = true
#     traverse(dihd.a3)

#     # return coordinates
#     state
# end

@inline function Base.max(state::State, comp::Symbol)
    x = getproperty(state, comp)
    m = zero(eltype(x))
    for i = 1:state.size
        f = x[i,1]^2 + x[i,2]^2 + x[i,3]^2
        if f > m
            m = f
        end
    end
    return sqrt(m)
end


function measure(dihd::Dihedral, xyz::Matrix{Float64})

    @nexprs 3 u -> v12_u = xyz[dihd.a1, u] - xyz[dihd.a0, u]
    @nexprs 3 u -> v23_u = xyz[dihd.a2, u] - xyz[dihd.a1, u]
    @nexprs 3 u -> v34_u = xyz[dihd.a3, u] - xyz[dihd.a2, u]
    @cross u n_u v12_u v23_u
    @cross u m_u v23_u v34_u
    @cross u o_u   n_u   m_u
    x = @dot(u, o_u, v23_u)/sqrt(@dot(u, v23_u, v23_u))
    y = @dot(u, n_u, m_u) 
    atan(x,y)
end

measure(dihd::Dihedral, state::State) = measure(dihd, state.coords)