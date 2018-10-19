using ProtoSyn
using LinearAlgebra
using Printf

# ------------------------------------------------------------------------------------------------------------
#                     STEEPEST DESCENT MINIMIZATION OF THE NATIVE STRUCTURE 

#1. LOAD STATE -----------------------------------------------------------------------------------------------
state = Common.load_from_pdb("data/1ctf_native.pdb")
state.energy = Forcefield.Amber.Energy()


#2. LOAD AMBER TOPOLOGY (BONDS, ANGLES, DIHEDRALS AND NON-BONDED) --------------------------------------------
topology = Forcefield.Amber.load_from_json("data/1ctf_amber_top.json")


#3. DEFINE THE EVALUATOR --------------------------------------------------------------------------------------
function my_evaluator!(st::Common.State, do_forces::Bool)
    energy = Forcefield.Amber.evaluate!(topology, st, cut_off=1.2, do_forces=do_forces)
    return energy
end


#4. DEFINE THE CALLBACKS --------------------------------------------------------------------------------------
callback1 = @Common.callback 100 function callback_status(step::Int64, st::Common.State, dr::Drivers.SteepestDescent.SteepestDescentDriver, args...)
    write(stdout, @sprintf "(SD) Step: %4d | Energy: %9.4f | Max Force: %9.4f | Gamma: %9.4f\n" step state.energy.eTotal args[1] args[2])
end
# callback1 = Common.CallbackObject(callback_status, 10)


file_xyz = open("out/trajectory.xyz", "w")
function callback_print(step::Int64, st::Common.State, dr::Drivers.SteepestDescent.SteepestDescentDriver, args...)
    Print.as_xyz(file_xyz, st, "Step: $step")
end
callback2 = Common.CallbackObject(1000, callback_print)


#5. DEFINE THE DRIVER PARAMETERS AND RUN THE SIMULATION -------------------------------------------------------
sd_driver = Drivers.SteepestDescent.SteepestDescentDriver(my_evaluator!, n_steps = 100000)
Drivers.SteepestDescent.run!(state, sd_driver, callback1, callback2)

close(file_xyz)