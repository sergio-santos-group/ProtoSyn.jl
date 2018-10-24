using ProtoSyn
using LinearAlgebra
using Printf

# ------------------------------------------------------------------------------------------------------------
#                                        TEST ENVIRONMENT 


file_xyz = open("out/trajectory.pdb", "w")
mc_topology = Aux.read_JSON("data/1ctf_mc_top.json")
dihedrals, residues = Common.load_topology(mc_topology)

index = 1
for file_i in ["1ctf.pdb", "1ctf.gro"]
    println("Reading file $file_i")
    string(file_i[end-2:end]) == "gro" ? state = Common.load_from_gro("data/$file_i") : state = Common.load_from_pdb("data/$file_i")
    state.energy = Forcefield.Amber.Energy()

    # Common.renumber_residues!(state.metadata.atoms)
    Common.fix_proline!(state, dihedrals)
    # Common.stretch_conformation!(state, dihedrals)
    # 4. DEFINE THE EVALUATOR -------------------------------------------------------------------------------------
    topology = Forcefield.Amber.load_from_json("data/1ctf_amber_top.json")
    function my_evaluator!(st::Common.State, do_forces::Bool)
        energy = Forcefield.Amber.evaluate!(topology, st, cut_off=1.2, do_forces=do_forces)
        return energy
    end


    # 5. DEFINE THE DRIVER -----------------------------------------------------------------------------------------
    sd_driver = Drivers.SteepestDescent.SteepestDescentDriver(my_evaluator!, n_steps = 1000)


    #4. DEFINE THE CALLBACK(S) -------------------------------------------------------------------------------------
    my_callback = @Common.callback 100 function callback_status(step::Int64, st::Common.State, dr::Drivers.SteepestDescent.SteepestDescentDriver, args...)
        write(stdout, @sprintf "(SD) Step: %4d | Energy: %9.4f | Max Force: %9.4f | Gamma: %9.4f\n" step state.energy.eTotal args[1] args[2])
        Print.as_pdb(file_xyz, st, step=step)
    end


    # 5. RUN THE SIMULATION -----------------------------------------------------------------------------------------
    Drivers.SteepestDescent.run!(state, sd_driver, my_callback)
    
    Print.as_pdb(file_xyz, state, step=index)
    exit(1)






    global index += 1
    # Print.as_pdb(file_xyz, state, step=index)
    # global index += 1
    Common.apply_ss!(state, dihedrals, "CEEEEEEECCCCHHHHHHHHHHHHCCCHHHHHHHHHCCCEEEEEEECHHHHHHHHHHHHHHCCEEEEC")
    Print.as_pdb(file_xyz, state, step=index)
    global index += 1
end
close(file_xyz)