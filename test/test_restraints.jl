using ProtoSyn
using Printf

function test_dihedralFBR_energy(workload::Int64 = 360)

    # Create 4 particles forming a dihedral
    state = Common.State(4)
    state.xyz = [0.1 0.1 0.1; 0.1 0.1 0.1; 0.1 0.1 0.1; 0.1 0.1 0.1]
    metadata = Common.Metadata(atoms = [Common.AtomMetadata(), Common.AtomMetadata(), Common.AtomMetadata(), Common.AtomMetadata()])

    # Create a makeshift Amber Forcefield
    init_angle = -180
    amber_ff = Forcefield.Amber.Topology()
    amber_ff.bonds = [Forcefield.Amber.HarmonicBond(1, 2, 1e5, 0.15), Forcefield.Amber.HarmonicBond(2, 3, 1e5, 0.15), Forcefield.Amber.HarmonicBond(3, 4, 1e5, 0.15)]
    amber_ff.angles = [Forcefield.Amber.HarmonicAngle(1, 2, 3, 120.0, deg2rad(120.0)), Forcefield.Amber.HarmonicAngle(2, 3, 4, 120.0, deg2rad(120.0))]
    amber_ff.dihedralsCos = [Forcefield.Amber.DihedralCos(1, 2, 3, 4, 40.0, deg2rad(init_angle), 2.0)]

    # Minimize the dihedral to a starting position, based on the Amber Forcefield
    function evaluate!(st::Common.State, do_forces::Bool)
        return Forcefield.Amber.evaluate!(amber_ff, st, cut_off=1.2, do_forces=do_forces)
    end
    sd_driver = Drivers.SteepestDescent.SteepestDescentDriver(evaluate!, n_steps = 100)
    Drivers.SteepestDescent.run!(state, sd_driver)
    
    # Define the DihedralFBR
    fbr = Forcefield.Restraints.DihedralFBR(1, 2, 3, 4, deg2rad(45), deg2rad(80), deg2rad(100), deg2rad(135), 1e3)

    Δangle = deg2rad(180/workload)
    for i in 1:workload
        Common.rotate_dihedral!(state.xyz, 2, 3, Δangle, Common.DIHEDRAL.omega, Int64[3, 4], nothing)
        println("$(init_angle + rad2deg(Δangle) * i) $(Forcefield.Restraints.evaluate!([fbr], state, do_forces=false))")
    end
end


# function test_dihedralFBR_forces()

#     # Create 4 particles forming a dihedral
#     state = Common.State(4)
#     state.xyz = [0.0 1.0 0.0; 0.0 0.0 0.0; 1.0 0.0 0.0; 1.0 1.0 0.0]
#     # metadata = Common.Metadata(atoms = [Common.AtomMetadata(), Common.AtomMetadata(), Common.AtomMetadata(), Common.AtomMetadata()])
    
#     # Define the DihedralFBR
#     fbr = Forcefield.Restraints.DihedralFBR(1, 2, 3, 4, deg2rad(45), deg2rad(80), deg2rad(100), deg2rad(135), 1e3)
    
#     Δangle = deg2rad(1e-5)
#     e1 = Forcefield.Restraints.evaluate!([fbr], state, do_forces=false)
#     Common.rotate_dihedral!(state.xyz, 2, 3, Δangle, Common.DIHEDRAL.omega, Int64[3, 4], nothing)
#     e2 = Forcefield.Restraints.evaluate!([fbr], state, do_forces=true)
#     f1 = abs(- (e2 - e1) / (Δangle))
#     f2 = abs(-sum(state.forces[1, :]))
#     Δf = abs(f1 - f2)
#     println(@sprintf "DIHEDRAL | F1(numeric) = %6.2e | F2(analytical) = %6.2e | ΔF = %5.2e | %s" f1 f2 Δf string(Δf < 0.0001))
# end


function test_distanceFBR_energy(workload::Int64 = 360)

    # Create 2 particles
    init_distance::Float64 = 10.0
    state = Common.State(2)
    state.xyz = [0.0 0.0 0.0; init_distance 0.0 0.0]
    metadata = Common.Metadata(atoms = [Common.AtomMetadata(), Common.AtomMetadata()])
    
    # Define the DihedralFBR
    fbr = Forcefield.Restraints.DistanceFBR(1, 2, 2.0, 4.0, 6.0, 8.0, 1e3)
    
    Δx::Float64 = init_distance/workload
    for i in 1:workload
        state.xyz[2, :] -= [Δx, 0.0, 0.0]
        println("$(init_distance - Δx * i) $(Forcefield.Restraints.evaluate!([fbr], state, do_forces=false))")
    end
end


function test_distanceFBR_forces()

    # Create 2 particles
    init_distance::Float64 = 4.0
    state = Common.State(2)
    state.xyz = [0.0 0.0 0.0; init_distance 0.0 0.0]
    metadata = Common.Metadata(atoms = [Common.AtomMetadata(), Common.AtomMetadata()])
    
    # Define the DihedralFBR
    fbr = Forcefield.Restraints.DistanceFBR(1, 2, 2.0, 4.0, 6.0, 8.0, 1e3)
    
    Δx::Float64 = 0.0000001
    e1 = Forcefield.Restraints.evaluate!([fbr], state, do_forces=false)
    state.xyz[2, :] -= [Δx, 0.0, 0.0]
    e2 = Forcefield.Restraints.evaluate!([fbr], state, do_forces=true)
    f1 = abs(- (e2 - e1) / (Δx))
    println(sum(state.forces[1, :], " ", sum(state.forces[2, :])))
    f2 = abs(-sum(state.forces[1, :]))
    Δf = abs(f1 - f2)
    println(@sprintf "DISTANCE | F1(numeric) = %10.8f | F2(analytical) = %10.8f | ΔF = %5.2e | %s" f1 f2 Δf string(Δf < 0.0001))
end


# test_dihedralFBR_energy()
# test_distanceFBR_energy()

function test_dihedralFBR_forces(workload::Int64 = 360)

    state = Common.State(4)
    state.xyz = [0.0 -1.0 0.0; 0.0 0.0 0.0; 1.0 0.0 0.0; 1.0 1.0 0.0]
    metadata = Common.Metadata(atoms = [Common.AtomMetadata(), Common.AtomMetadata(), Common.AtomMetadata(), Common.AtomMetadata()])
    fbr = Forcefield.Restraints.DihedralFBR(1, 2, 3, 4, deg2rad(-90), deg2rad(-45), deg2rad(45), deg2rad(90), 1e3)
    
    for θ in [-180.0, 180.0]
        θ = deg2rad(θ)
        # Rotate atom 1, measure new energy and forces.
        state.xyz[1, 2] = cos(θ)
        state.xyz[1, 3] = sin(θ)
        fill!(state.forces, zero(Float64))
        e1 = Forcefield.Restraints.evaluate!([fbr], state, do_forces=true)
        f1 = deepcopy(state.forces)
        f2 = zeros(size(f1, 1), size(f1, 2))

        ϵ = 1e-8
        for atom_index in 1:size(state.xyz, 1)
            # Forces need to be measured and compared by moving one component at a time
            for component in [1, 2, 3]
                backup_component = state.xyz[atom_index, component]
                state.xyz[atom_index, component] += ϵ
                e2 = Forcefield.Restraints.evaluate!([fbr], state, do_forces=false)
                f2[atom_index, component] = - (e2 - e1) / ϵ
                state.xyz[atom_index, component] = backup_component
            end
        end
        # Print.as_pdb(stdout, state, metadata)
        println("$(floor(Int64, rad2deg(θ))) $e1 $(maximum(f2-f1)) $f1 $f2")
        # println("$(floor(Int64, rad2deg(θ))) $e1 $(maximum(f2-f1))")
    end

end

test_dihedralFBR_forces()
# test_distanceFBR_forces()