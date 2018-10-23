using ProtoSyn
using LinearAlgebra
using Printf

# -------------------------------------------------------------------------------------------------------------
#                     DIHEDRAL MOVEMENTS ON THE BACKBONE WITH BLOCKED SECONDARY STRUCTURES 

# 1. LOAD STATE -----------------------------------------------------------------------------------------------
state = Common.load_from_pdb("data/1ctf.pdb")
state.energy = Forcefield.Amber.Energy()


# 2. APPLY SECONDARY STRUCTURE --------------------------------------------------------------------------------
mc_topology = Aux.read_JSON("data/1ctf_mc_top.json")
dihedrals, residues = Common.load_topology(mc_topology)
bb_dihedrals = filter(x -> x.dtype < Common.omega, dihedrals)
Common.apply_initial_conf!(state, bb_dihedrals)


# 3. FIX PROLINES ---------------------------------------------------------------------------------------------
Common.fix_proline!(state, dihedrals);


# 4. DEFINE THE SAMPLER ----------------------------------------------------------------------------------------
bb_nb_dihedrals = filter(x -> x.dtype < Common.omega && Int(x.residue.ss) < 1, dihedrals)
dihedral_mutator = Mutators.Dihedral.DihedralMutator(bb_nb_dihedrals, randn, 0.1, 1.0)
function my_sampler(st::Common.State)
    Mutators.Dihedral.run!(st, dihedral_mutator)
end


# 5. DEFINE THE EVALUATOR --------------------------------------------------------------------------------------
topology = Forcefield.Amber.load_from_json("data/1ctf_amber_top.json")
function my_evaluator!(st::Common.State, do_forces::Bool)
    energy = Forcefield.Amber.evaluate!(topology, st, cut_off=1.2, do_forces=do_forces)
    return energy
end


# 6. DEFINE THE DRIVER -----------------------------------------------------------------------------------------
mc_driver = Drivers.MonteCarlo.MonteCarloDriver(my_sampler, my_evaluator!, temperature=1.0, n_steps=1000)


# 7. DEFINE THE CALLBACKS --------------------------------------------------------------------------------------
callback1 = @Common.callback 10 function callback_status(step::Int64, st::Common.State, dr::Drivers.MonteCarlo.MonteCarloDriver, args...)
    write(stdout, @sprintf "(MC) Step: %4d | Energy: %9.4f\n" step state.energy.eTotal)
end

file_xyz = open("out/trajectory.xyz", "w")
function callback_print(step::Int64, st::Common.State, dr::Drivers.MonteCarlo.MonteCarloDriver, args...)
    Print.as_xyz(file_xyz, st, "Step: $step")
end
callback2 = Common.CallbackObject(100, callback_print)


# 8. RUN THE SIMULATION ---------------------------------------------------------------------------------------
Drivers.MonteCarlo.run!(state, mc_driver, callback1, callback2)
close(file_xyz)