using ProtoSyn
using LinearAlgebra
using Printf

# -------------------------------------------------------------------------------------------------------------
#   DIHEDRAL+CRANKSHAFT MOVEMENTS ON THE BACKBONE WITH BLOCKED SECONDARY STRUCTURES AND STEP SIZE ADJUSTMENT

# 1. LOAD STATE -----------------------------------------------------------------------------------------------
state = Common.load_from_pdb("data/1ctf.pdb")
state.energy = Forcefield.Amber.Energy()


# 2. APPLY SECONDARY STRUCTURE --------------------------------------------------------------------------------
mc_topology = Aux.read_JSON("data/1ctf_mc_top.json")
dihedrals, residues = Common.load_topology(mc_topology)
bb_dihedrals = filter(x -> x.dtype < Common.omega, dihedrals)
Common.apply_ss!(state, bb_dihedrals, "CEEEEEEECCCCHHHHHHHHHHHHCCCHHHHHHHHHCCCEEEEEEECHHHHHHHHHHHHHHCCEEEEC")


# 3. FIX PROLINES ---------------------------------------------------------------------------------------------
Common.fix_proline!(state, dihedrals);


# 4. DEFINE THE SAMPLER ----------------------------------------------------------------------------------------
bb_nb_dihedrals = filter(x -> x.dtype < Common.omega && Int(x.residue.ss) < 1, dihedrals)
dihedral_mutator = Mutators.Dihedral.DihedralMutator(bb_nb_dihedrals, () -> (rand() * 2 - 1 * dihedral_mutator.step_size), 0.01, 0.01)
crankshaft_mutator = Mutators.Crankshaft.CrankshaftMutator(bb_nb_dihedrals, () -> (rand() * 2 - 1 * crankshaft_mutator.step_size), 0.005, 0.01)
function my_sampler!(st::Common.State)
    Mutators.Dihedral.run!(st, dihedral_mutator)
    Mutators.Crankshaft.run!(st, crankshaft_mutator)
end


# 5. DEFINE THE EVALUATOR --------------------------------------------------------------------------------------
topology = Forcefield.Amber.load_from_json("data/1ctf_amber_top.json")
function my_evaluator!(st::Common.State, do_forces::Bool)
    energy = Forcefield.Amber.evaluate!(topology, st, cut_off=1.2, do_forces=do_forces)
    return energy
end


# 6. DEFINE THE DRIVER -----------------------------------------------------------------------------------------
mc_driver = Drivers.MonteCarlo.MonteCarloDriver(my_sampler!, my_evaluator!, temperature=500.0, n_steps=1000)


# 7. DEFINE THE CALLBACKS --------------------------------------------------------------------------------------
callback1 = @Common.callback 100 function cb_status(step::Int64, st::Common.State, dr::Drivers.MonteCarlo.MonteCarloDriver, args...)
    write(stdout, @sprintf "(MC) Step: %4d | Energy: %.4ef | D step: %.3e | CS step: %.3e | AR: %.2f\n" step state.energy.eTotal dihedral_mutator.step_size crankshaft_mutator.step_size args[1])
end

file_xyz = open("out/trajectory.pdb", "w")
callback2 = @Common.callback 100 function cb_print(step::Int64, st::Common.State, dr::Drivers.MonteCarlo.MonteCarloDriver, args...)
    Print.as_pdb(file_xyz, st, step = step)
end

callback3 = @Common.callback 1 function cb_adjust_step_size(step::Int64, st::Common.State, dr::Drivers.MonteCarlo.MonteCarloDriver, args...)
    acceptance_ratio = 0.2
    if args[1] > acceptance_ratio
        dihedral_mutator.step_size * 1.05 < π ? dihedral_mutator.step_size *= 1.05 : dihedral_mutator.step_size = π
        crankshaft_mutator.step_size * 1.1 < π ? crankshaft_mutator.step_size *= 1.1 : crankshaft_mutator.step_size = π
    else
        dihedral_mutator.step_size * 0.95 > 1e-5 ? dihedral_mutator.step_size *= 0.95 : dihedral_mutator.step_size = 1e-5
        crankshaft_mutator.step_size * 0.9 > 1e-5 ? crankshaft_mutator.step_size *= 0.9 : crankshaft_mutator.step_size = 1e-5
    end
end


# 8. RUN THE SIMULATION ---------------------------------------------------------------------------------------
Drivers.MonteCarlo.run!(state, mc_driver, callback1, callback2, callback3)
close(file_xyz)