using Base.Cartesian
using ...ProtoSyn: @dot, @cross

#Calculators.eval!(s::State, c::ComponentContainer{T}, f::Bool=true) where

function Calculators.eval!(state::State, bonds::ComponentContainer{HarmonicBond}, do_forces::Bool=false)
    
    energy = 0.0
    forces = state.forces
    coords = state.coords
    
    @inbounds for bond in bonds.items
        
        a1 = bond.a1
        a2 = bond.a2
        type = bond.type
        
        @nexprs 3 u -> v12_u = coords[a2,u] - coords[a1,u]
        
        d12 = sqrt(@dot(u, v12_u, v12_u))
            
        dr = d12 - type.r
        energy += type.k * dr^2
        
        if do_forces
            fc = type.k * dr / d12
            @nexprs 3 u -> forces[a1, u] += fc * v12_u
            @nexprs 3 u -> forces[a2, u] -= fc * v12_u
        end

    end
    state.energy.HarmonicBond = energy * 0.5
    state.energy.total = energy * 0.5
end



function Calculators.eval!(state::State, angles::ComponentContainer{HarmonicAngle}, do_forces::Bool=false)
    
    energy::Float64 = 0.0
    forces = state.forces
    coords = state.coords

    @inbounds for angle in angles.items

        a1 = angle.a1
        a2 = angle.a2
        a3 = angle.a3
        type = angle.type

        @nexprs 3 u -> v21_u = coords[a1,u] - coords[a2,u]
        @nexprs 3 u -> v23_u = coords[a3,u] - coords[a2,u]

        # 1.0/(v21⋅v21)
        # inv_d21sq  = 1.0 / @ncall(3, +, u -> v21_u^2)
        inv_d21sq  = 1.0 / @dot(u, v21_u, v21_u)
        
        # 1.0/(v23⋅v23)
        # inv_d23sq  = 1.0 / @ncall(3, +, u -> v23_u^2)
        inv_d23sq  = 1.0 / @dot(u, v23_u, v23_u)
        
        # cos(θ) = v21⋅v23/(|v21|*|v23|)
        # ctheta = @ncall(3, +, u -> v21_u * v23_u)
        ctheta = @dot(u, v21_u, v23_u)

        inv_d21xd23 = sqrt(inv_d21sq * inv_d23sq)
        ctheta *= inv_d21xd23
        
        dtheta  = acos(ctheta) - type.θ
        energy += type.k * dtheta^2
        
        if do_forces
            # prevent division by zero
            fc = max(1.0 - ctheta^2, 1e-8)
            fc = dtheta * type.k / sqrt(fc)
            # fc = dtheta * type.k / sqrt(1.0 - ctheta*ctheta)
            # println("Δθ = $dtheta, k=$(type.k), cos(θ)=$ctheta, fc=$fc, E=$energy")
            
            @nexprs 3 u -> begin
                f1 = fc*(ctheta * v21_u*inv_d21sq - v23_u*inv_d21xd23)
                f3 = fc*(ctheta * v23_u*inv_d23sq - v21_u*inv_d21xd23)
                forces[a1, u] -= f1
                forces[a2, u] += (f1 + f3)
                forces[a3, u] -= f3
            end
            
        end
    end
    
    state.energy.HarmonicAngle = energy * 0.5
    state.energy.total = energy * 0.5
end




function Calculators.eval!(state::State, dihedrals::ComponentContainer{CosineDihedral}, do_forces::Bool=false)
    
    energy = 0.0
    forces = state.forces
    coords = state.coords
    
    @inbounds for dihd in dihedrals.items

        a1 = dihd.a1
        a2 = dihd.a2
        a3 = dihd.a3
        a4 = dihd.a4
        type = dihd.type

        @nexprs 3 u -> v12_u = coords[a2,u] - coords[a1,u]
        @nexprs 3 u -> v32_u = coords[a2,u] - coords[a3,u]
        @nexprs 3 u -> v34_u = coords[a4,u] - coords[a3,u]
        # d32sq = @ncall 3 (+) u -> v32_u^2
        d32sq = @dot u v32_u v32_u

        # vm = v12 × v32 = cross(v12, v32)
        # vm_1 = v12_2*v32_3 - v12_3*v32_2
        # vm_2 = v12_3*v32_1 - v12_1*v32_3
        # vm_3 = v12_1*v32_2 - v12_2*v32_1
        @cross u vm_u v12_u v32_u
        
        
        # vn = v32 × v34 = cross(v32, v34)
        # vn_1 = v32_2*v34_3 - v32_3*v34_2
        # vn_2 = v32_3*v34_1 - v32_1*v34_3
        # vn_3 = v32_1*v34_2 - v32_2*v34_1
        @cross u vn_u v32_u v34_u

        d32 = sqrt(d32sq)
        # v12⋅vn
        # d12n = @ncall 3 (+) u -> v12_u * vn_u
        d12n = @dot u v12_u vn_u
        # vm⋅vn
        # dmn  = @ncall 3 (+) u ->  vm_u * vn_u
        dmn  = @dot u vm_u vn_u

        phi = atan(d32 * d12n, dmn)
        energy += type.k * (1.0 + cos(type.n * phi - type.θ))
        
        if do_forces
            
            dVdphi_x_d32 = type.k * type.n * sin(type.θ - type.n * phi) * d32
            
            # v34⋅v32
            # d3432 = @ncall 3 (+) u -> v34_u * v32_u
            d3432 = @dot(u, v34_u, v32_u)
            # v12⋅v32
            # d1232 = @ncall 3 (+) u -> v12_u * v32_u
            d1232 = @dot(u, v12_u, v32_u)
            # vm⋅vm
            # dmm = @ncall 3 (+) u -> vm_u^2
            dmm = @dot(u, vm_u, vm_u)
            # vn⋅vn
            # dnn = @ncall 3 (+) u -> vn_u^2
            dnn = @dot(u, vn_u, vn_u)
            
            c3_1 = (d3432/d32sq - 1.0)
            c3_2 = (d1232/d32sq)

            c1   = -dVdphi_x_d32 / dmm
            c4   =  dVdphi_x_d32 / dnn

            @nexprs 3 u -> begin
                f1 = vm_u * c1
                f4 = vn_u * c4
                f3 = f4 * c3_1 - f1 * c3_2
                forces[a1, u] -= f1
                forces[a2, u] -= (-f1 - f3 - f4)
                forces[a3, u] -= f3
                forces[a4, u] -= f4
            end
        end

    end
    
    state.energy.CosineDihedral = energy
end


function Calculators.eval!(state::State, dihedrals::ComponentContainer{NCosineDihedral}, do_forces::Bool=false)
    
    energy = 0.0
    forces = state.forces
    coords = state.coords
    
    #@inbounds
    for dihd in dihedrals.items

        a1 = dihd.a1
        a2 = dihd.a2
        a3 = dihd.a3
        a4 = dihd.a4
        type = dihd.type
        
        @nexprs 3 u -> v12_u = coords[a2,u] - coords[a1,u]
        @nexprs 3 u -> v32_u = coords[a2,u] - coords[a3,u]
        @nexprs 3 u -> v34_u = coords[a4,u] - coords[a3,u]
        d32sq = @dot u v32_u v32_u

        # vm = v12 × v32 = cross(v12, v32)
        vm_1 = v12_2*v32_3 - v12_3*v32_2
        vm_2 = v12_3*v32_1 - v12_1*v32_3
        vm_3 = v12_1*v32_2 - v12_2*v32_1
        
        
        # vn = v32 × v34 = cross(v32, v34)
        vn_1 = v32_2*v34_3 - v32_3*v34_2
        vn_2 = v32_3*v34_1 - v32_1*v34_3
        vn_3 = v32_1*v34_2 - v32_2*v34_1

        d32 = sqrt(d32sq)
        # dot(v12, vn)
        #d12n = @ncall 3 (+) u -> v12_u * vn_u
        d12n = @dot u v12_u vn_u
        # dot(vm, vn)
        # dmn  = @ncall 3 (+) u ->  vm_u * vn_u
        dmn  = @dot u vm_u vn_u

        phi = atan(d32 * d12n, dmn)
        

        # for i=1:dihd.size
        for i in eachindex(type.k)
            energy += type.k[i] * (1.0 + cos(type.n[i] * phi - type.θ[i]))
        end

        if do_forces
            
            # dot(v34, v32)
            # d3432 = @ncall 3 (+) u -> v34_u * v32_u
            d3432 = @dot u v34_u v32_u
            # dot(v12, v32)
            # d1232 = @ncall 3 (+) u -> v12_u * v32_u
            d1232 = @dot u v12_u v32_u
            # dot(vm, vm)
            # dmm = @ncall 3 (+) u -> vm_u^2
            dmm = @dot u vm_u vm_u
            # dot(vn, vn)
            # dnn = @ncall 3 (+) u -> vn_u^2
            dnn = @dot u vn_u vn_u
            
            c3_1 = (d3432/d32sq - 1.0)
            c3_2 = (d1232/d32sq)
            
            for i in eachindex(type.k)
                dVdphi_x_d32 = type.k[i] * type.n[i] * sin(type.θ[i] - type.n[i] * phi) * d32
                c1   = -dVdphi_x_d32 / dmm
                c4   =  dVdphi_x_d32 / dnn

                @nexprs 3 u -> begin
                    f1 = vm_u * c1
                    f4 = vn_u * c4
                    f3 = f4 * c3_1 - f1 * c3_2
                    forces[a1, u] -= f1
                    forces[a2, u] -= (-f1 - f3 - f4)
                    forces[a3, u] -= f3
                    forces[a4, u] -= f4
                end
            end

        end

    end
    
    state.energy.NCosineDihedral = energy
end



# function evaluate!(state::State, angles::Vector{HarmonicAngle}, do_forces::Bool=false)
#     energy = 0.0
#     forces = state.forces
#     coords = state.coords
    
#     v21 = [0.0,0.0,0.0]
#     v23 = [0.0,0.0,0.0]
    
#     @inbounds for angle in angles

#         a1 = angle.a1
#         a2 = angle.a2
#         a3 = angle.a3

#         ctheta = 0.0
#         d21Sq = 0.0
#         d23Sq = 0.0

#         for i=1:3
#             v21[i] = coords[a1,i] - coords[a2,i]
#             v23[i] = coords[a3,i] - coords[a2,i]
#             d21Sq += v21[i]^2
#             d23Sq += v23[i]^2
#             ctheta += v21[i]*v23[i]
#         end

#         d21xd23 = sqrt(d21Sq * d23Sq)
#         ctheta /= d21xd23
#         dtheta  = acos(ctheta) - angle.θ

#         energy += angle.k * dtheta^2
        

#         if do_forces
#             fc = dtheta * angle.k / sqrt(1.0 - ctheta * ctheta)
#             for i=1:3
#                 f1 = fc*(ctheta * v21[i]/d21Sq - v23[i]/d21xd23)
#                 f3 = fc*(ctheta * v23[i]/d23Sq - v21[i]/d21xd23)
                
#                 forces[a1, i] += f1
#                 forces[a3, i] += f3
#                 forces[a2, i] -= (f1 + f3)
#             end

#         end
#     end
    
#     energy *= 0.5

#     energy
# end




# function ProtoSyn.eval!(state::State, top::Topology, do_forces::Bool=false)
#     for field in fieldnames(Topology)
#         ProtoSyn.eval!(state, getfield(top, field))
#     end
# end