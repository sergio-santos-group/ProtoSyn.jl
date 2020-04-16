using Base.Cartesian
using LinearAlgebra: norm, cross

export rotmat, rotmat!, rotate!, measure

@inline function rotmat!(R::Matrix{T}, axis::Vector{T}, θ::T) where {T<:AbstractFloat}
    # q0 = cos(0.5 * θ)
    # q1, q2, q3 = sin(0.5 * θ) * axis ./ norm(axis)
    λ  = sin(θ/2) / norm(axis)
    q0 = cos(θ/2)
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

@inline rotmat(axis::Vector{T}, θ::T) where {T<:AbstractFloat} = begin
    rotmat!(zeros(T, 3, 3), axis, θ)
end

rotate!(state::State{T}, mol::Molecule, block::AbstractAxisRotatableBlock, θ::T) where {T<:AbstractFloat} = begin
    update!(mol.traverse_state, block, mol.bonds)
    
    # axis of rotation along vector a1->a2
    axis = T[0, 0, 0]
    xyz = state.coords
    @nexprs 3 k -> axis[k] = xyz[k, block.a2] - xyz[k, block.a1]
    
    # get rotation matrix
    rmat = T[0 0 0;0 0 0; 0 0 0]
    rotmat!(rmat, axis, θ)
    # rmat = rotmat(axis, θ)
    
    # perform rotation and return rotated state
    rotate!(state.coords, rmat, block.a1, mol.traverse_state)
    return state
end


rotate!(xyz::Matrix{T}, rmat::Matrix{T}, pivot::Int, ts::TraverseState) where {T<:AbstractFloat} = begin
    indices = ts.indices
    @inbounds begin
        @nexprs 3 k -> o_k = xyz[k, pivot]
        for j = 1:ts.size
            i = indices[j]
            @nexprs 3 k -> x_k = xyz[k, i] - o_k
            @nexprs 3 k -> xyz[k, i] = rmat[k,1]*x_1 + rmat[k,2]*x_2 + rmat[k,3]*x_3 + o_k
            # @nexprs 3 k -> x_k = xyz[i,k] - o_k
            # @nexprs 3 k -> xyz[i,k] = rmat[k,1]*x_1 + rmat[k,2]*x_2 + rmat[k,3]*x_3 + o_k
        end
    end
    xyz
end

rotate!(xyz::Matrix{T}, rmat::Matrix{T}) where {T<:AbstractFloat} = begin
    @inbounds begin
        for i = 1:size(xyz, 2)
            @nexprs 3 k -> x_k = xyz[k, i]
            @nexprs 3 k -> xyz[k, i] = rmat[k,1]*x_1 + rmat[k,2]*x_2 + rmat[k,3]*x_3
        end
    end
    xyz
end

function measure(dihd::Dihedral, xyz::Matrix{Float64})

    @nexprs 3 u -> v12_u = xyz[u, dihd.a1] - xyz[u, dihd.a0]
    @nexprs 3 u -> v23_u = xyz[u, dihd.a2] - xyz[u, dihd.a1]
    @nexprs 3 u -> v34_u = xyz[u, dihd.a3] - xyz[u, dihd.a2]
    @cross u n_u v12_u v23_u
    @cross u m_u v23_u v34_u
    @cross u o_u   n_u   m_u
    x = @dot(u, o_u, v23_u)/sqrt(@dot(u, v23_u, v23_u))
    y = @dot(u, n_u, m_u) 
    atan(x,y)
end

measure(dihd::Dihedral, state::State) = measure(dihd, state.coords)

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
