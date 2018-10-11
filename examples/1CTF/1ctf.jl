using ProtoSyn

file_xyz = open("teste_output.xyz", "w")

state = Common.load_from_pdb("mol2.pdb")
state.energy = Forcefield.Energy()

topology = Forcefield.load_from_json("1ctf.json")

mc_topology = Aux.read_JSON("1ctf_mc_top.json")
dihedrals, residues = Mutators.Dihedral.load_topology(mc_topology)

#Add dummy dihedral for first C_alpha
insert!(dihedrals, 1, Mutators.Dihedral.NewDihedral(-1, 2, 3, length(residues[1].atoms) - 1, residues[1].atoms[4:end], residues[1], "PHI"))

#Mask only the PHI dihedrals
phi_dihedrals = filter(x -> x.dtype == "PHI", dihedrals)

steepest_descent_params = Drivers.SteepestDescent.ConfigParameters(n_steps = 300, log_freq = 100)
diehdral_params = Mutators.Dihedral.ConfigParameters(p_mut = 5e-3)
crankshaft_params = Mutators.Crankshaft.ConfigParameters(p_mut = 1e-3)
montecarlo_params = Drivers.MonteCarlo.ConfigParameters(n_steps = 100000, log_freq = 1, step_size = 0.33)

function my_evaluator1!(st::Common.State, do_forces::Bool)
    return Forcefield.evalenergy!(topology, st, cut_off=1.2, do_forces=do_forces)
end

function my_evaluator2!(st::Common.State, do_forces::Bool)
    energy::Float64 = 0.0
    fill!(state.forces, 0.0)
    energy += Forcefield.evaluate!(topology.bonds, st, do_forces=do_forces)
    energy += Forcefield.evaluate!(topology.angles, st, do_forces=do_forces)
    energy += Forcefield.evaluate!(topology.dihedralsCos, st, do_forces=do_forces)
    return energy
end

function my_sampler!(st::Common.State)
    diehdral_params.step_size = montecarlo_params.step_size
    crankshaft_params.step_size = montecarlo_params.step_size
    Mutators.Dihedral.run!(st, dihedrals, diehdral_params, () -> randn() * π * diehdral_params.step_size)
    Mutators.Crankshaft.run!(st, phi_dihedrals, crankshaft_params, () -> randn() * π * crankshaft_params.step_size)
    Drivers.SteepestDescent.run!(st, my_evaluator2!, steepest_descent_params)
end

function my_callback(st::Common.State, step::Int)
    Print.as_xyz(st, ostream = file_xyz, title = "Step $step")
end

my_callback(state, 0)
@time Drivers.MonteCarlo.run!(state, my_sampler!, my_evaluator1!, montecarlo_params, callback = my_callback)
close(file_xyz)