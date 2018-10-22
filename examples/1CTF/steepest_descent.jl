using ProtoSyn
using LinearAlgebra
using Printf

# ------------------------------------------------------------------------------------------------------------
#                     STEEPEST DESCENT MINIMIZATION OF THE NATIVE STRUCTURE 

#1. LOAD STATE -----------------------------------------------------------------------------------------------
file_xyz = open("out/trajectory.pdb", "w")
state = Common.load_from_pdb("data/1ctf.pdb")
state.energy = Forcefield.Amber.Energy()

#2. Apply initial state ?
mc_topology = Aux.read_JSON("data/1ctf_mc_top.json")
dihedrals, residues = Common.load_topology(mc_topology)
bb_dihedrals = filter(x -> x.dtype < Common.omega, dihedrals)
Common.apply_dihedrals_from_file(state, bb_dihedrals, "data/1ctf_native.pdb")
Print.as_pdb(file_xyz, state, step=2)

#2.1 FIX PROLINE !
state = Common.load_from_pdb("data/1ctf.pdb")
state.energy = Forcefield.Amber.Energy()
Print.as_pdb(file_xyz, state, step=3)
# for dihedral in dihedrals
#     if dihedral.residue.name == "P" && dihedral.dtype == Common.phi
#         r::Common.Residue = dihedral.residue
#         Common.rotate_dihedral!(state.xyz, 565, 567, deg2rad(180), Common.phi, r.atoms, r)
#         insert!(dihedral.movable, 1, 569)
#         insert!(dihedral.movable, 1, 568)
#     end

# end
Common.fix_proline(state, dihedrals)
Print.as_pdb(file_xyz, state, step=4)
bb_dihedrals = filter(x -> x.dtype < Common.omega, dihedrals)
Common.apply_dihedrals_from_file(state, bb_dihedrals, "1ctf.pdb")
Print.as_pdb(file_xyz, state, step=5)


#2. LOAD AMBER TOPOLOGY (BONDS, ANGLES, DIHEDRALS AND NON-BONDED) --------------------------------------------
topology = Forcefield.Amber.load_from_json("data/1ctf_amber_top.json")


#3. DEFINE THE EVALUATOR --------------------------------------------------------------------------------------
function my_evaluator!(st::Common.State, do_forces::Bool)
    energy = Forcefield.Amber.evaluate!(topology, st, cut_off=1.2, do_forces=do_forces)
    return energy
end


#4. DEFINE THE CALLBACK(S) --------------------------------------------------------------------------------------
my_callback = @Common.callback 100 function callback_status(step::Int64, st::Common.State, dr::Drivers.SteepestDescent.SteepestDescentDriver, args...)
    write(stdout, @sprintf "(SD) Step: %4d | Energy: %9.4f | Max Force: %9.4f | Gamma: %9.4f\n" step state.energy.eTotal args[1] args[2])
    Print.as_pdb(file_xyz, st, step=step)
end


#5. DEFINE THE DRIVER PARAMETERS AND RUN THE SIMULATION -------------------------------------------------------
sd_driver = Drivers.SteepestDescent.SteepestDescentDriver(my_evaluator!, n_steps = 1000)
Drivers.SteepestDescent.run!(state, sd_driver, my_callback)

close(file_xyz)