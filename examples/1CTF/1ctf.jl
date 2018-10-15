using ProtoSyn
using LinearAlgebra

file_xyz = open("teste_output.pdb", "w")
state = Common.load_from_pdb("mol.pdb")

state.energy = Forcefield.Energy()

topology = Forcefield.load_from_json("1ctf.json")

mc_topology = Aux.read_JSON("1ctf_mc_top.json")
dihedrals, residues = Mutators.Dihedral.load_topology(mc_topology)

#Add dummy dihedral for first C_alpha
insert!(dihedrals, 1, Mutators.Dihedral.NewDihedral(-1, 2, 3, length(residues[1].atoms) - 1, residues[1].atoms[4:end], residues[1], "PHI"))

#Mask only the PHI dihedrals
phi_dihedrals = filter(x -> x.dtype == "PHI", dihedrals)
bb_dihedrals = filter(x -> x.dtype in ["PHI", "PSI"], dihedrals)

steepest_descent_params = Drivers.SteepestDescent.ConfigParameters(n_steps = 100, log_freq = 100, callback_freq = 50)
diehdral_params = Mutators.Dihedral.ConfigParameters(p_mut = 5e-5)
crankshaft_params = Mutators.Crankshaft.ConfigParameters(p_mut = 5e-4)
montecarlo_params = Drivers.MonteCarlo.ConfigParameters(n_steps = 1000000, log_freq = 100, step_size = 1.0, temperature = 150.0)

# ---
struct SolvPair
    i::Int64
    coef::Float64
end

d = Dict("Q"=>-3.5,"W"=>-0.9,"T"=>-0.7,"C"=>2.5,"P"=>-1.6,"V"=>4.2,"L"=>3.8,"M"=>1.9,"N"=>-3.5,"H"=>-3.2,"A"=>1.8,"D"=>-3.5,"G"=>-0.4,"E"=>-3.5,"Y"=>-1.3,"I"=>4.5,"S"=>-0.8,"K"=>-3.9,"R"=>-4.5,"F"=>2.8)
tmp = map(x -> SolvPair(x.a3, d[string(x.residue.name[1])]), phi_dihedrals)

function calc_eSol(st::Common.State)
    n_atoms = length(tmp)
    l = zeros(Float64, 3)
    e_sol::Float64 = 0.0
    sum_f::Float64 = 0.0
    for i in 1:(n_atoms - 1)
        Ai = @view st.xyz[tmp[i].i,:]
        sum_f = 0.0
        for j in (i+1):n_atoms
            Aj = @view st.xyz[tmp[j].i,:]
            @. l[:] = Aj - Ai
            dIJ = norm(l)
            f = 1.0 / (1.0 + exp(-(12.0-dIJ) / 0.4))
            sum_f += f
        end
        if ((sum_f < 21.0) && (tmp[i].coef > 0.0)) || ((sum_f > 21.0) && (tmp[i].coef < 0.0))
            e_sol += tmp[i].coef * (21.0 - sum_f)
        end
    end
    return e_sol
end

# ---

function my_evaluator1!(st::Common.State, do_forces::Bool)
    energy = Forcefield.evaluate!(topology, st, cut_off=1.2, do_forces=do_forces)
    eSol = calc_eSol(st)
    # println("eSol: $eSol")
    return energy + eSol
    # return Forcefield.evaluate!(topology, st, cut_off=1.2, do_forces=do_forces)
end

function my_evaluator2!(st::Common.State, do_forces::Bool)
    energy::Float64 = 0.0
    fill!(state.forces, 0.0)
    energy += Forcefield.evaluate!(topology.bonds, st, do_forces=do_forces)
    energy += Forcefield.evaluate!(topology.angles, st, do_forces=do_forces)
    energy += Forcefield.evaluate!(topology.dihedralsCos, st, do_forces=do_forces)
    return energy
end

movements = open("mov.txt", "w")
function my_sampler!(st::Common.State)
    # diehdral_params.step_size = montecarlo_params.step_size
    # crankshaft_params.step_size = montecarlo_params.step_size
    Mutators.Dihedral.run!(st, bb_dihedrals, diehdral_params, () -> randn() * π * diehdral_params.step_size)
    Mutators.Crankshaft.run!(st, phi_dihedrals, crankshaft_params, () -> randn() * π * crankshaft_params.step_size)
    Drivers.SteepestDescent.run!(st, my_evaluator2!, steepest_descent_params)
end

function my_callback(st::Common.State, step::Int)
    Print.as_pdb(st, ostream = file_xyz, title = "Step $step")
end

# Drivers.SteepestDescent.run!(state, my_evaluator1!, steepest_descent_params)
my_callback(state, 0)
@time Drivers.MonteCarlo.run!(state, my_sampler!, my_evaluator1!, montecarlo_params, callback = my_callback)
close(file_xyz)
close(movements)