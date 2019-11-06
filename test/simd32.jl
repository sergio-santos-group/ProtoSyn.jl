using Base.Cartesian
using SIMD



mutable struct VerletList
    capacity::Int               # max size of list 
    cutoff::Float32             # cutoff (nm)
    buffer::Float32             # buffer (nm)
    offset::Vector{Int}    # pointer particle->list position
    list::Vector{Int}      # nblist
    # _
end
VerletList(n::Int) = VerletList(n, -1.0, 0.0, zeros(Int, n), zeros(Int, n))


function update!(vlist::VerletList, coords::Matrix{Float32})
    if vlist.cutoff < 0.0
        return vlist
    end
    
    ptr = 1
    natoms = size(coords, 1)
    cutsq = (vlist.cutoff + vlist.buffer)^2

    for i = 1:natoms
        vlist.offset[i] = ptr
        lastptr = ptr
        for j = (i+1):natoms
            @nexprs 3 u -> vij_u = coords[i,u] - coords[j,u]
            dijsq = vij_1^2 + vij_2^2 + vij_3^2
            if dijsq < cutsq
                vlist.list[ptr] = j
                ptr += 1
                if ptr == vlist.capacity
                    resize!(vlist, vlist.capacity + 1000)
                end
            end
        end

        # pad section
        nfills = ptr - lastptr
        if nfills % 4 > 0
            nf = 4 - nfills%4
            #last = vlist.list[ptr-1]
            for i=1:nf
                vlist.list[ptr] = natoms+1#last
                ptr += 1
                if ptr == vlist.capacity
                    resize!(vlist, vlist.capacity + 1000)
                end
            end
        end
        
        
        vlist.list[ptr] = -1
        ptr += 1

        if (i < natoms) && (ptr == vlist.capacity)
            resize!(vlist, vlist.capacity + 1000)
        end

    end
    vlist
end
Base.resize!(vlist::VerletList, n::Int) = begin
    resize!(vlist.list, n)
    vlist.capacity = n
    vlist
end

#-------------------------


# function lj(xyz::Matrix{T}) where {T}
#     n = size(xyz, 1)
#     energy = zero(T)
#     σ = one(T)

#     @inbounds for i=1:n-1
#         @nexprs 3 u 
#         for j = i+1:n
#             @nexprs 3 u -> dij_u = xyz[j,u] - xyz[i,u]
#             dij_sq = dij_1^2 + dij_2^2 + dij_3^2
#             lj2 = σ/dij_sq
#             lj6 = lj2*lj2*lj2
#             energy += lj6*lj6 - lj6
#         end
#     end
#     energy
# end

function lj(xyz::Matrix{T}) where {T}
    n = size(xyz, 1)
    energy = zero(T)
    σ = one(T)

    @inbounds for i=1:n-1
        @nexprs 3 u -> xi_u = xyz[i,u]
        
        for j = i+1:n
            # @nexprs 3 u -> dij_u = xyz[j,u] - xyz[i,u]
            @nexprs 3 u -> dij_u = xyz[j,u] - xi_u

            dij_sq = dij_1^2 + dij_2^2 + dij_3^2
            lj2 = σ/dij_sq
            lj6 = lj2*lj2*lj2
            energy += lj6*lj6 - lj6
        end
    end
    energy
end


function lj_swp(xyz::Matrix{T}) where {T}
    n = size(xyz, 1)
    energy = zero(T)
    σ = one(T)

    @inbounds for i=1:n-1

        @nexprs 3 u -> xi_u = xyz[i,u]
        @nexprs 3 u -> dij_u = xyz[i+1,u] - xi_u
        dij_sq = dij_1^2 + dij_2^2 + dij_3^2

        if i+2 < n
            for j = i+2:n
                lj2 = σ/dij_sq
                lj6 = lj2*lj2*lj2
                energy += lj6*lj6 - lj6
                
                @nexprs 3 u -> dij_u = xyz[j,u] - xi_u
                dij_sq = dij_1^2 + dij_2^2 + dij_3^2
            end
        end
        lj2 = σ/dij_sq
        lj6 = lj2*lj2*lj2
        energy += lj6*lj6 - lj6
    end
    energy
end



function lj(xyz::Matrix{T}, vlist::VerletList) where {T}
    n = size(xyz, 1)
    energy = zero(T)
    σ = one(T)

    @inbounds for i=1:n-1
        ptr =vlist.offset[i]
        while vlist.list[ptr] > 0
            j = vlist.list[ptr]
            ptr += 1

            @nexprs 3 u -> dij_u = xyz[j,u] - xyz[i,u]
            dij_sq = dij_1^2 + dij_2^2 + dij_3^2

            lj2 = σ/dij_sq
            lj6 = lj2*lj2*lj2
            energy += lj6*lj6 - lj6
        end
    end
    energy
end



function lj_simd(xyz::Matrix{T}, ϵ::Vector{T},::Type{Vec{N,T}}) where {N, T}
    n = size(xyz, 1)
    
    x = view(xyz, :, 1)
    y = view(xyz, :, 2)
    z = view(xyz, :, 3)
    
    σ = one(T)
    energy = zero(T)

    @inbounds for i=1:n-1
        xi = Vec{N,T}(x[i])
        yi = Vec{N,T}(y[i])
        zi = Vec{N,T}(z[i])
        for j=i+1:4:n
            xj = vload(Vec{N,T}, x, j)
            yj = vload(Vec{N,T}, y, j)
            zj = vload(Vec{N,T}, z, j)

            xij_sq = (xj - xi)^2
            yij_sq = (yj - yi)^2
            zij_sq = (zj - zi)^2
            dij_sq = xij_sq + yij_sq + zij_sq
            
            lj2 = σ/dij_sq
            lj6 = lj2*lj2*lj2
            m = vload(Vec{N,T}, ϵ, j)
            energy += sum(m*(lj6*lj6 - lj6))
        end
    end
    energy
end



function lj_simd2(xyz::Matrix{T}, vlist::VerletList, ϵ::Vector{T}, xt::Vector{T}, ::Type{Vec{N,T}}) where {N, T}
    n = size(xyz, 1)
    m = 2*n
    o = 3*n
    σ =one(T)
    energy = zero(T)
    
    x = view(xyz, :, 1)
    y = view(xyz, :, 2)
    z = view(xyz, :, 3)
    @inbounds for i=1:n-1
        
        ptr = vlist.offset[i]
        if vlist.list[ptr] < 1
            continue
        end
        
        xi = Vec{N,T}(xyz[i,1])
        yi = Vec{N,T}(xyz[i,2])
        zi = Vec{N,T}(xyz[i,3])
        

        # j = vlist.list[ptr]
        # k = 1
        # while j > 0
        #     # xt[  k] = xyz[j, 1]
        #     # xt[n+k] = xyz[j, 2]
        #     # xt[m+k] = xyz[j, 3]
        #     # xt[o+k] = ϵ[j]
        #     ptr += 1
        #     j = vlist.list[ptr]
        #     k += 1
        # end

        #for j=1:N:k-1
        #println("$i  ->  ptr=$(ptr), list[ptr]=$(vlist.list[ptr])")
        while vlist.list[ptr] > 0
            #j = vlist.list[ptr]
        #    ptr += 4

            idx = vload(Vec{N,Int}, vlist.list, ptr)
            ptr += N
            xj = vgather(x,   idx)
            yj = vgather(y, n+idx)
            zj = vgather(z, m+idx)
            mj = vgather(ϵ,   idx)

            # xj = vload(Vec{N,T}, xt,   j)
            # yj = vload(Vec{N,T}, xt, n+j)
            # zj = vload(Vec{N,T}, xt, m+j)
            # mj = vload(Vec{N,T}, xt, o+j)

            xij_sq = (xj - xi)^2
            yij_sq = (yj - yi)^2
            zij_sq = (zj - zi)^2
            dij_sq = xij_sq + yij_sq + zij_sq
            # vij_sq = (kj - ki)^2
            # dij_sq = sum(vij_sq)
            # println("$i  $j   ", dij_sq)
            
            lj2 = σ/dij_sq
            lj6 = lj2*lj2*lj2
            #m = vload(V464, ϵ, j)
            # mask = vgather(ϵ, idx)
            energy += sum(mj*(lj6*lj6 - lj6))
            # energy += sum(lj6*lj6 - lj6)
        end
        # println(ptr)
    end
    energy
end




# function lj_simd3(xyz::Vector{T}, vlist::VerletList, ϵ::Vector{T}, ::Type{Vec{N,T}}) where {N, T}
function lj_simd3(n::Int, xyz::Vector{T}, vlist::VerletList, ::Type{Vec{N,T}}) where {N, T}
    
    σ = one(T)
    energy = zero(T)
    m1 = (0, 4, 1, 5)
    m2 = (2, 6, 3, 7)
    mx = (0, 1, 4, 5)
    my = (2, 3, 6, 7)

    # n = NParticles
    # println(vlist)
    @inbounds for i=1:n-1
        
        ptr = vlist.offset[i]
        if vlist.list[ptr] < 1
            continue
        end
        stop = vlist.offset[i+1]-2

        vi = vload(Vec{N,T}, xyz, (i-1)*4+1)
        
        for j = ptr:N:stop
            
            jj = (vlist.list[  j]-1)*4+1
            v1 = vload(Vec{N,T}, xyz, jj)

            jj = (vlist.list[j+1]-1)*4+1
            v2 = vload(Vec{N,T}, xyz, jj)
            
            jj = (vlist.list[j+2]-1)*4+1
            v3 = vload(Vec{N,T}, xyz, jj)
            
            jj = (vlist.list[j+3]-1)*4+1
            v4 = vload(Vec{N,T}, xyz, jj)
            
            
            vi1 = v1 - vi # xi1, yi1, zi1, wi1
            vi2 = v2 - vi # xi2, yi2, zi2, wi2
            vi3 = v3 - vi # xi3, yi3, zi3, wi3
            vi4 = v4 - vi # xi4, yi4, zi4, wi4

            vs12 = shufflevector(vi1, vi2, Val{m1})    # xi1, xi2, yi1, yi2
            vs34 = shufflevector(vi3, vi4, Val{m1})    # xi3, xi4, yi3, yi4
            
            xij = shufflevector(vs12, vs34, Val{mx})   # xi1, xi2, xi3, xi4
            yij = shufflevector(vs12, vs34, Val{my})   # yi1, yi2, yi3, yi4
            

            vs12 = shufflevector(vi1, vi2, Val{m2})    # zi1, zi2, ?i1, ?i2
            vs34 = shufflevector(vi3, vi4, Val{m2})    # zi3, zi4, ?i3, ?i4

            zij = shufflevector(vs12, vs34, Val{mx})   # zi1, zi2, zi3, zi4

            dij_sq = xij^2 + yij^2 + zij^2

            lj2 = σ/dij_sq
            lj6 = lj2*lj2*lj2
            #m = vload(V464, ϵ, j)
            # mask = vgather(ϵ, idx)
            energy += sum(lj6*lj6 - lj6)
        end
    end
    energy
end





NParticles = 10_000
xyz = 4*rand(Float32, NParticles,3)


vl = VerletList(NParticles)
vl.cutoff = 1.2
vl.buffer = 0.2

update!(vl, xyz)

function work_lj(n)
    e = zero(eltype(xyz))
    for i=1:n
        e = lj(xyz)
    end
    println("  > ", e)
end

function work_lj_swp(n)
    e = zero(eltype(xyz))
    for i=1:n
        e = lj_swp(xyz)
    end
    println("  > ", e)
end

function work_lj_vl(n)
    e = zero(eltype(xyz))
    for i=1:n
        e = lj(xyz, vl)
    end
    println("  > ", e)
end

function work_lj_simd(n)
    e = zero(eltype(xyz))
    
    ϵ = ones(Float32, NParticles+8)
    ϵ[NParticles+1:end] .= zero(eltype(xyz))


    for i=1:n
        e = lj_simd(xyz, ϵ, Vec{4,Float32})
    end
    println("  > ", e)
end


function work_lj_vl_simd(n)
    e = zero(eltype(xyz))

    for i=1:n
        xyz2 = vcat(xyz', zeros(Float32, 1, NParticles))[:]
        e = lj_simd3(NParticles, xyz2, vl, Vec{8,Float32})
    end
    println("  > ", e)

end


Nreps = 20
println("SERIAL")
@time work_lj(1)
@time work_lj(Nreps)

println("SERIAL SWP")
@time work_lj_swp(1)
@time work_lj_swp(Nreps)

println("SERIAL VLIST")
@time work_lj_vl(1)
@time work_lj_vl(Nreps)

println("SIMD")
@time work_lj_simd(1)
@time work_lj_simd(Nreps)

println("SIMD VLIST")
@time work_lj_vl_simd(1)
@time work_lj_vl_simd(Nreps)

