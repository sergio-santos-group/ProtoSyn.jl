using Base.Cartesian
using SIMD



mutable struct VerletList
    capacity::Int               # max size of list 
    cutoff::Float64             # cutoff (nm)
    buffer::Float64             # buffer (nm)
    offset::Vector{Int}    # pointer particle->list position
    list::Vector{Int}      # nblist
    # _
end
VerletList(n::Int) = VerletList(n, -1.0, 0.0, zeros(Int, n), zeros(Int, n))


function update!(vlist::VerletList, coords::Matrix{Float64})
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


function lj(xyz::Matrix{Float64})
    n = size(xyz, 1)
    energy = 0.0
    @inbounds for i=1:n-1
        for j = i+1:n
            @nexprs 3 u -> dij_u = xyz[j,u] - xyz[i,u]
            dij_sq = dij_1^2 + dij_2^2 + dij_3^2
            # println("$i  $j   ", dij_sq)
            #if dij_sq > 100.0
            #    continue
            #end

            lj2 = 1.0/dij_sq
            lj6 = lj2*lj2*lj2
            energy += lj6*lj6 - lj6
        end
    end
    energy
end


function lj(xyz::Matrix{Float64}, vlist::VerletList, ::Matrix{Float64})
    n = size(xyz, 1)
    energy = 0.0
    @inbounds for i=1:n-1
        ptr =vlist.offset[i]
        while vlist.list[ptr] > 0
            j = vlist.list[ptr]
            ptr += 1

            @nexprs 3 u -> dij_u = xyz[j,u] - xyz[i,u]
            dij_sq = dij_1^2 + dij_2^2 + dij_3^2
            # println("$i  $j   ", dij_sq)
            #if dij_sq > 100.0
            #    continue
            #end

            lj2 = 1.0/dij_sq
            lj6 = lj2*lj2*lj2
            energy += lj6*lj6 - lj6
        end
    end
    energy
end




const V464 = Vec{4,Float64}
const V264 = Vec{2,Float64}
V864 = Vec{8,Float64}

function lj_simd(xyz::Matrix{Float64})
    n = size(xyz, 1)
    
    x = view(xyz, :, 1)
    y = view(xyz, :, 2)
    z = view(xyz, :, 3)
    
    ϵ = ones(n+4)
    ϵ[n+1:end] .= 0.0

    energy = 0.0
    @inbounds for i=1:n-1
        xi = V464(x[i])
        yi = V464(y[i])
        zi = V464(z[i])
        for j=i+1:4:n
            xj = vload(V464, x, j)
            yj = vload(V464, y, j)
            zj = vload(V464, z, j)

            xij_sq = (xj - xi)^2
            yij_sq = (yj - yi)^2
            zij_sq = (zj - zi)^2
            dij_sq = xij_sq + yij_sq + zij_sq
            # println("$i  $j   ", dij_sq)
            
            lj2 = 1.0/dij_sq
            lj6 = lj2*lj2*lj2
            m = vload(V464, ϵ, j)
            energy += sum(m*(lj6*lj6 - lj6))
        end
    end
    energy
end



# function lj_simd(xyz::Matrix{Float64}, vlist::VerletList, tmp::Matrix{Float64})
#     n = size(xyz, 1)
#     #println("aqui")
    

#     # tmp = similar(xyz)

#     x  = view(xyz, :, 1)
#     y  = view(xyz, :, 2)
#     z  = view(xyz, :, 3)

#     xt = view(tmp, :, 1)
#     yt = view(tmp, :, 2)
#     zt = view(tmp, :, 3)

#     energy = 0.0
#     @inbounds for i=1:n-1
#         xi = V464(x[i])
#         yi = V464(y[i])
#         zi = V464(z[i])

#         ptr = vlist.offset[i]
#         if vlist.list[ptr] < 1
#             continue
#         end
#         stop = vlist.offset[i+1]-2

#         #copyto!(tmp, view(xyz, idx, :))

#         nj = stop-ptr+1
#         #vtmp = view(tmp, 1:nj, :)
#         v = view(vlist.list, ptr:stop)
#         copyto!(xt, view(x, v))
#         copyto!(yt, view(y, v))
#         copyto!(zt, view(z, v))
#         # copyto!(vtmp, view(xyz, v, :))
        
#         #println(ptr, ",", stop, ",",nj, ",")
#         for j = 1:4:nj
#             xj = vload(V464, xt, j)
#             yj = vload(V464, yt, j)
#             zj = vload(V464, zt, j)

#             xij_sq = (xj - xi)^2
#             yij_sq = (yj - yi)^2
#             zij_sq = (zj - zi)^2
#             dij_sq = xij_sq + yij_sq + zij_sq
#             # println("$i  $j   ", dij_sq)
            
#             lj2 = 1.0/dij_sq
#             lj6 = lj2*lj2*lj2
#             energy += sum(lj6*lj6 - lj6)
#             # m = vload(V464, ϵ, j)
#             # energy += sum(m*(lj6*lj6 - lj6))
#         end

#     end
#     energy
# end

# function lj_simd(x::Vector{Float64}, y::Vector{Float64}, z::Vector{Float64}, vlist::VerletList)
function lj_simd2(xyz::Matrix{Float64}, vlist::VerletList, ϵ::Vector{Float64}, xt::Vector{Float64})
    n = size(xyz, 1)
    m = 2*n

    # ϵ = ones(n+4)
    # ϵ[n+1:end] .= 0.0
    # xyz = view(xyz, :)
    # println("ola")
    energy = 0.0
    @inbounds for i=1:n-1
        xi = V464(xyz[i,1])
        yi = V464(xyz[i,2])
        zi = V464(xyz[i,3])
        # xi = V464(xyz[  i])
        # yi = V464(xyz[n+i])
        # zi = V464(xyz[m+i])
        
        ptr = vlist.offset[i]
        if vlist.list[ptr] < 1
            continue
        end

        j = vlist.list[ptr]
        k = 1
        while j > 0
            xt[  k] = xyz[j, 1]
            xt[n+k] = xyz[j, 2]
            xt[m+k] = xyz[j, 3]
            ptr += 1
            j = vlist.list[ptr]
            k += 1
        end

        for j=1:4:k-1
        #println("$i  ->  ptr=$(ptr), list[ptr]=$(vlist.list[ptr])")
        #while vlist.list[ptr] > 0
        #    j = vlist.list[ptr]
        #    ptr += 4

            # idx = vload(Vec{4,Int}, vlist.list, ptr-4)
            
            # xj = vgather(xyz,   idx)
            # yj = vgather(xyz, n+idx)
            # zj = vgather(xyz, m+idx)

            xj = vload(V464, xt,   j)
            yj = vload(V464, xt, n+j)
            zj = vload(V464, xt, m+j)

            xij_sq = (xj - xi)^2
            yij_sq = (yj - yi)^2
            zij_sq = (zj - zi)^2
            dij_sq = xij_sq + yij_sq + zij_sq
            # vij_sq = (kj - ki)^2
            # dij_sq = sum(vij_sq)
            # println("$i  $j   ", dij_sq)
            
            lj2 = 1.0/dij_sq
            lj6 = lj2*lj2*lj2
            #m = vload(V464, ϵ, j)
            # mask = vgather(ϵ, idx)
            # energy += sum(mask*(lj6*lj6 - lj6))
            energy += sum(lj6*lj6 - lj6)
        end
    end
    energy
end

# function lj_simd2(xyz::Matrix{T}, ::Type{Vec{N,T}}) where {N,T}
    
#     x = view(xyz, :, 1)
#     y = view(xyz, :, 2)
#     z = view(xyz, :, 3)
    
#     n = size(xyz, 1)
#     ϵ = ones(T, n+4)
    
#     ϵ[n+1:end] .= zero(T)

#     energy = zero(T)
#     sigma = one(T)

#     @inbounds for i = 1:n-1
#         xi = Vec{N,T}(x[i])
#         yi = Vec{N,T}(y[i])
#         zi = Vec{N,T}(z[i])
        
#         for j = i+1:N:n
#             xj = vload(Vec{N,T}, x, j)
#             yj = vload(Vec{N,T}, y, j)
#             zj = vload(Vec{N,T}, z, j)

#             xij_sq = (xj - xi)^2
#             yij_sq = (yj - yi)^2
#             zij_sq = (zj - zi)^2
#             dij_sq = xij_sq + yij_sq + zij_sq
#             # println("$i  $j   ", dij_sq)
            
#             lj2 = sigma/dij_sq
#             lj6 = lj2*lj2*lj2

#             m = vload(Vec{N,T}, ϵ, j)
#             energy += sum(m*(lj6*lj6 - lj6))

#             # fc = 24.0 * m * (lj6 - 2.0 * lj6 * lj6) / d12sq

#         end
#     end
#     energy
# end

NParticles = 10_000
xyz = 4*rand(NParticles,3)

# x = xyz[:,1]
# y = xyz[:,2]
# z = xyz[:,3]

# xyz2 = vcat(xyz', zeros(1, NParticles))[:]

# xyz = [ 1.20246   4.94875  0.05359
#         2.72789   2.19132  4.79887
#         0.119921  3.72672  1.66502
#         1.1002    2.08762  0.565046
#         0.225788  4.51229  0.7467]

 vl = VerletList(NParticles)
 vl.cutoff = 1.2
 vl.buffer = 0.2

 update!(vl, xyz)


function work(n::Int, f::F) where {F <: Function}
    e = 0.0
    for i=1:n
        e = f(xyz)
    end
    println("  > ", e)
end

function work(n::Int, f::F, l::VerletList) where {F <: Function}
    tmp = similar(xyz)
    e = 0.0
    for i=1:n
        e = f(xyz, l, tmp)
    end
    println("  > ", e)
end

function w2(n)
    e = 0.0
    N = size(xyz,1)
    ϵ = ones(N+4)
    ϵ[n+1:end] .= 0.0
    xt = zeros(3*N+3)
    for i=1:n
        e = lj_simd2(xyz, vl, ϵ, xt)
        # e = lj_simd(x, y, z, vl)
    end
    println("  > ", e)
end

Nreps = 10
println("SERIAL")
@time work(1, lj)
@time work(Nreps, lj)

println("SERIAL VLIST")
@time work(1, lj, vl)
@time work(Nreps, lj, vl)

println("SIMD")
@time work(1, lj_simd)
@time work(Nreps, lj_simd)

println("SIMD VLIST")
@time w2(1)
@time w2(Nreps)
# @time work(1, lj_simd, vl)
# @time work(Nreps, lj_simd, vl)


# # println(lj(xyz))