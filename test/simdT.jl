using Base.Cartesian
using SIMD

using Random
Random.seed!(1234)

using BenchmarkTools



macro reduce(n::Int, op::Symbol, ex::Expr)
    esc(:(@ncall($n, $op, $ex)))
end




mutable struct VerletList
    size::Int              # number of particles the list refers to 
    capacity::Int          # capacity of list 
    cutoff::Float64        # cutoff (nm)
    buffer::Float64        # buffer (nm)
    offset::Vector{Int}    # pointer particle->list position
    list::Vector{Int}      # nblist
end

VerletList(n::Int) = VerletList(n, n, -1.0, 0.0, zeros(Int, n), zeros(Int, n))


function update!(vlist::VerletList, coords::Matrix{T}) where {T<:AbstractFloat}
    
    @assert vlist.size == max(size(coords)...) "incompatible sizes"
    # println(size(coords))

    # effective cutoff squared
    reff = vlist.cutoff + vlist.buffer
    cutsq = convert(T, reff*reff)
    
    offset = 1
    natoms = vlist.size
    @inbounds for i = 1:natoms

        vlist.offset[i] = offset
        # @nexprs 3 u -> xi_u = coords[i,u]
        @nexprs 3 u -> xi_u = coords[u,i]

        for j = (i+1):natoms

            @nexprs 3 u -> vij_u = coords[u,j] - xi_u
            # @nexprs 3 u -> vij_u = coords[j,u] - xi_u
            dij_sq = @reduce 3 (+) u -> vij_u*vij_u
            
            if dij_sq < cutsq
                vlist.list[offset] = j
                offset += 1
                if offset == vlist.capacity
                    resize!(vlist, vlist.capacity + natoms)
                end
            end
        end

        vlist.list[offset] = -1
        offset += 1
        if (i < natoms) && (offset == vlist.capacity)
            resize!(vlist, vlist.capacity + natoms)
        end

    end
    
    vlist
end

Base.resize!(vlist::VerletList, n::Int) = begin
    resize!(vlist.list, n)
    vlist.capacity = n
    vlist
end

# ==============================================================================
struct State{T<:AbstractFloat}
    size::Int
    coords::Matrix{T}
    forces::Matrix{T}
end


# ==============================================================================
include("simd_serial_aos.jl")
include("simd_serial_soa.jl")
include("simd_soa.jl")
include("simd_aos.jl")
include("simd_coords.jl")

# ==============================================================================

# # ==============================================================================
# function seria2l(xyz::Matrix{T}, forces::Matrix{T}, vlist::VerletList) where {T <: AbstractFloat}
    
#     n = max(size(xyz)...)
#     energy = zero(T)
#     σ = one(T)
#     cutsq = convert(T, vlist.cutoff * vlist.cutoff)
    
#     @inbounds for i=1:n-1
        
#         @nexprs 3 u -> xi_u = xyz[i,u]
#         @nexprs 3 u -> fi_u = zero(T)
        
#         ptr = vlist.offset[i]
#         while vlist.list[ptr] > 0
#             j = vlist.list[ptr]
#             ptr += 1

#             @nexprs 3 u -> rij_u = xyz[j,u] - xi_u
        
#             dij_sq = rij_1^2 + rij_2^2 + rij_3^2
#             if dij_sq > cutsq
#                 continue
#             end
            
#             lj2 = σ/dij_sq
#             lj6 = lj2*lj2*lj2
#             energy += lj6*lj6 - lj6

#             # do forces
#             fc = T(24.0) * (lj6 - T(2.0) * lj6 * lj6) / dij_sq
#             @nexprs 3 u -> begin
#                 fi_u += fc * rij_u
#                 forces[j, u] -= fc * rij_u
#             end

#         end
#         @nexprs 3 u -> forces[i, u] += fi_u
#     end
    
#     energy
# end


# # ==============================================================================
# function simd(xyz::Matrix{T}, vlist::VerletList{T}, ::Type{Vec{4,T}}) where {T}

#     N = 4

#     σ = one(T)
#     n = size(xyz, 1)
#     energy = zero(T)
    
#     xyz = vcat(xyz', zeros(T, 1, n))[:]
    
#     @inbounds for i=1:n-1

#         ptr = vlist.offset[i]
#         if vlist.list[ptr] < 1
#             continue
#         end
#         stop = vlist.offset[i+1]-2

#         # load point i
#         vi = vload(Vec{N,T}, xyz, 4(i-1)+1)


#         #println(i)
        
#         xi = Vec{N,T}(x[i])
#         yi = Vec{N,T}(y[i])
#         zi = Vec{N,T}(z[i])

#         # number of i-j interactions
#         ii = i+1
#         ni = n - ii + 1

#         # if number of interactions is greater than the
#         # SIMD vector size, we can perform, at least, one cycle of SIMD
#         # vectorized calculations
#         padlen = rem(ni,N)
#         if ni >= N
#             for j = ii:N:n-padlen
#                 #println("  s:",j:N+j-1)
#                 xj = vload(Vec{N,T}, x, j)
#                 yj = vload(Vec{N,T}, y, j)
#                 zj = vload(Vec{N,T}, z, j)
            
#                 xij_sq = (xj - xi)^2
#                 yij_sq = (yj - yi)^2
#                 zij_sq = (zj - zi)^2
#                 dij_sq = xij_sq + yij_sq + zij_sq

#                 lj2 = σ/dij_sq
#                 lj6 = lj2*lj2*lj2

#                 energy += sum(lj6*lj6 - lj6)
#             end
#         end
#         for j = n-padlen+1:n
#             #println("  i:",j)
#             @nexprs 3 u -> sdij_u = xyz[j,u] - xyz[i,u]
#             sdij_sq = sdij_1^2 + sdij_2^2 + sdij_3^2
#             slj2 = σ / sdij_sq
#             slj6 = slj2 * slj2 * slj2
#             energy += slj6*slj6 - slj6
#         end
#     end
#     energy
# end

# # ==============================================================================
# function simd(xyz::Matrix{T}, forces::Matrix{T}, ::Type{Vec{N,T}}) where {N, T}

#     σ = one(T)
#     n = size(xyz, 1)
#     energy = zero(T)
    
#     x = view(xyz, :, 1)
#     y = view(xyz, :, 2)
#     z = view(xyz, :, 3)
    
#     fx = view(forces, :, 1)
#     fy = view(forces, :, 2)
#     fz = view(forces, :, 3)
#     lane = VecRange{N}(0)

#     @inbounds for i=1:n-1
#         #println(i)
        
#         xi = Vec{N,T}(x[i])
#         yi = Vec{N,T}(y[i])
#         zi = Vec{N,T}(z[i])

#         fxi = Vec{N,T}(0) # fx[i])
#         fyi = Vec{N,T}(0) # fy[i])
#         fzi = Vec{N,T}(0) # fz[i])

#         # number of i-j interactions
#         ii = i+1
#         ni = n - ii + 1

#         # if number of interactions is greater than the
#         # SIMD vector size, we can perform, at least, one cycle of SIMD
#         # vectorized calculations
#         padlen = rem(ni,N)
#         if ni >= N
#             for j = ii:N:n-padlen
#                 #println("  s:",j:N+j-1)
                
#                 xj = vload(Vec{N,T}, x, j)
#                 yj = vload(Vec{N,T}, y, j)
#                 zj = vload(Vec{N,T}, z, j)
            
#                 rxij = xj - xi
#                 ryij = yj - yi
#                 rzij = zj - zi
                
#                 dij_sq = rxij^2 + ryij^2 + rzij^2

#                 lj2 = σ/dij_sq
#                 lj6 = lj2*lj2*lj2

#                 energy += sum(lj6*lj6 - lj6)

#                 # do forces
#                 fc = T(24.0) * (lj6 - T(2.0) * lj6 * lj6) / dij_sq

#                 fxi += fc * rxij
#                 fyi += fc * ryij
#                 fzi += fc * rzij
                
#                 lanej = lane+j
#                 fx[lanej] -= fc * rxij
#                 fy[lanej] -= fc * ryij
#                 fz[lanej] -= fc * rzij

#             end
#         end
#         for j = n-padlen+1:n
#             #println("  i:",j)
#             @nexprs 3 u -> srij_u = xyz[j,u] - xyz[i,u]
#             sdij_sq = srij_1^2 + srij_2^2 + srij_3^2
#             slj2 = σ / sdij_sq
#             slj6 = slj2 * slj2 * slj2
#             energy += slj6*slj6 - slj6

#             sfc = T(24.0) * (slj6 - T(2.0) * slj6 * slj6) / sdij_sq
#             @nexprs 3 u -> begin
#                 forces[i, u] += sfc * srij_u
#                 forces[j, u] -= sfc * srij_u
#             end
#         end

#         fx[i] += sum(fxi)
#         fy[i] += sum(fyi)
#         fz[i] += sum(fzi)

#     end
#     energy
# end

# ==============================================================================


# https://discourse.julialang.org/t/how-to-code-faster-parallel-for-loop/19760

# ==============================================================================
function simd(xyz::Matrix{T}, forces::Matrix{T}, vlist::VerletList, ::Type{Vec{4,T}}) where {T}
    
    n = size(xyz, 1)
    σ = one(T)
    energy = zero(T)
    
    # xyz = vcat(xyz', zeros(T, 1, n))[:]
    xyz = view(vcat(xyz', zeros(T, 1, n)), :)
    forces = zeros(T, 4n)
    lane = VecRange{4}(0)

    m1 = (0, 4, 1, 5)
    m2 = (2, 6, 3, 7)
    mx = (0, 1, 4, 5)
    my = (2, 3, 6, 7)

    @inbounds for i=1:n-1

        ptr = vlist.offset[i]
        if vlist.list[ptr] < 1
            continue
        end
        stop = vlist.offset[i+1]-2

        fi = Vec{4,T}(0)

        vi = vload(Vec{4,T}, xyz, 4(i-1)+1)

        for j = ptr:4:stop-4

            j1 = 4(vlist.list[j  ]-1)+1
            j2 = 4(vlist.list[j+1]-1)+1
            j3 = 4(vlist.list[j+2]-1)+1
            j4 = 4(vlist.list[j+3]-1)+1

            v1 = vload(Vec{4,T}, xyz, j1)
            v2 = vload(Vec{4,T}, xyz, j2)
            v3 = vload(Vec{4,T}, xyz, j3)
            v4 = vload(Vec{4,T}, xyz, j4)
            
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



            # do forces
            fc = T(24.0) * (lj6 - T(2.0) * lj6 * lj6) / dij_sq
            fi += fc[1] * vi1
            fi += fc[2] * vi2
            fi += fc[3] * vi3
            fi += fc[4] * vi4
            forces[lane+j1] -= fc[1] * vi1
            forces[lane+j2] -= fc[2] * vi2
            forces[lane+j3] -= fc[3] * vi3
            forces[lane+j4] -= fc[4] * vi4
        end

        forces[lane +  4(i-1)+1] += fi
    end
    energy
end


# # ==============================================================================
# function simd_4N(xyz::Matrix{T}, ::Type{Vec{4,T}}) where {T}
    
#     N = 4
#     σ = one(T)
#     n = size(xyz, 2)
#     energy = zero(T)
    
#     xyz = view(xyz,:)

#     @inbounds for i=1:n-1
#         #println(i)
        
#         xi = vload(Vec{N,T}, xyz, 4(i-1)+1)

#         for j=i+1:n
#             xj = vload(Vec{N,T}, xyz, 4(j-1)+1)
#             vij = xj - xi
#             dij_sq = sum(vij^2)
            
#             lj2 = σ/dij_sq
#             lj6 = lj2*lj2*lj2
#             energy += lj6*lj6 - lj6
#         end
#     end
#     energy
# end

# ==============================================================================
const FTYPE = Float32
# const NATOMS = 50_000
#const L = convert(FTYPE, 10.0)
x = generate_cubic_lattice(10.0, 31, LATTICE.face_centered, false)
NATOMS = size(x, 2)
println("NATOMS = $NATOMS")

xyz = convert(Matrix{FTYPE}, x)
forces = zeros(FTYPE, 4, NATOMS)
# xyz = L*rand(FTYPE, 4, NATOMS)
# forces = zeros(FTYPE, 4, NATOMS)
# xyz[4,:] .= FTYPE(0)


state1 = State(NATOMS, xyz, forces)
state2 = State(NATOMS, Array(xyz'), Array(forces'))
# xyz = collect(xyz')

list = VerletList(NATOMS)
list.cutoff = 3.0
list.buffer = 0.0

@time update!(list, xyz)
# @time update!(list, xyz)

println("SERIAL AOS w/ FORCES")
#@time
serial_aos(state1, list.cutoff, Val{true})
@time serial_aos(state1, list.cutoff, Val{true})

println("SERIAL AOS+VLIST w/ FORCES")
#@time
serial_aos(state1, list, Val{true})
@time serial_aos(state1, list, Val{true})

# println("SERIAL AOS w/o FORCES")
# @time  serial_aos(state1, list.cutoff, Val{false})
# @time serial_aos(state1, list.cutoff, Val{false})

println("SERIAL SOA w/ FORCES")
#@time
serial_soa(state2, list.cutoff, Val{true})
@time serial_soa(state2, list.cutoff, Val{true})

println("SERIAL SOA+VLIST w/ FORCES")
#@time
serial_soa(state2, list, Val{true})
@time serial_soa(state2, list, Val{true})

# println("SERIAL SOA w/o FORCES")
# @time  serial_soa(state2, list.cutoff, Val{false})
# @time serial_soa(state2, list.cutoff, Val{false})

println("SIMD-4 SOA w/ FORCES")
#@time
simd_soa(state2, list.cutoff, Val{true}, Vec{4,FTYPE})
@time simd_soa(state2, list.cutoff, Val{true}, Vec{4,FTYPE})

println("SIMD-8 SOA w/ FORCES")
#@time
simd_soa(state2, list.cutoff, Val{true}, Vec{8,FTYPE})
@time simd_soa(state2, list.cutoff, Val{true}, Vec{8,FTYPE})

# println("SIMD SOA w/o FORCES")
# @time  simd_soa(state2, list.cutoff, Val{false}, Vec{4,FTYPE})
# @btime simd_soa(state2, list.cutoff, Val{false}, Vec{4,FTYPE})

# println(state1)

println("SIMD-4 AOS+VLIST w/ FORCES")
#@time
simd_aos(state1, list, Val{true}, Vec{4,FTYPE})
@time simd_aos(state1, list, Val{true}, Vec{4,FTYPE})

println("SIMD-8 AOS+VLIST w/ FORCES")
#@time
simd_aos(state1, list, Val{true}, Vec{8,FTYPE})
@time simd_aos(state1, list, Val{true}, Vec{8,FTYPE})

# println(state1)
# println("SERIAL w/ LIST")
# @time  serial(xyz, forces, list)
# @btime serial(xyz, forces, list)

# println("SIMD")
# @time  simd(xyz, forces, Vec{8, FTYPE})
# @btime simd(xyz, forces, Vec{8, FTYPE})

# println("SIMD w/ LIST")
# @time  simd(xyz, forces, list, Vec{4, FTYPE})
# @btime simd(xyz, forces, list, Vec{4, FTYPE})

