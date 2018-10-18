using ProtoSyn
using LinearAlgebra
using Printf

# ------------------------------------------------------------------------------------

file_xyz = open("teste_output.xyz", "w")

#1. LOAD STATE
state = Common.load_from_pdb("mol.pdb")
state.energy = Forcefield.Amber.Energy()

#2. LOAD TOPOLOGY (DIHEDRALS AND RESIDUES)
mc_topology = Aux.read_JSON("1ctf_mc_top.json")
dihedrals, residues = Common.load_topology(mc_topology)

#3. LOAD AMBER TOPOLOGY (BONDS, ANGLES, DIHEDRALS AND NON-BONDED)
topology = Forcefield.Amber.load_from_json("1ctf.json")

#4. FILTER DIHEDRALS
bb_dihedrals = filter(x -> x.dtype < Common.omega, dihedrals)
bb_nb_dihedrals = filter(x -> x.dtype < Common.omega && Int(x.residue.ss) < 1, dihedrals)

#5. APPLY INITIAL CONFORMATION
Common.apply_initial_conf!(state, bb_dihedrals)

#6. DEFINE THE SAMPLER
dihedral_mutator = Mutators.Dihedral.DihedralMutator(bb_nb_dihedrals, 0.1, randn, 1.0)
function my_sampler(st::Common.State)
    Mutators.Dihedral.run!(st, dihedral_mutator)
end

#7. DEFINE THE EVALUATOR
function my_evaluator!(st::Common.State, do_forces::Bool)
    energy = Forcefield.Amber.evaluate!(topology, st, cut_off=1.2, do_forces=do_forces)
    return energy
end

#8. DEFINE THE CALLBACK
function callback(step::Int64, st::Common.State, dr::Drivers.MonteCarlo.MonteCarloDriver, args...)
    write(stdout, @sprintf "(MC) Step: %4d | Energy: %9.4f\n" step state.energy.eTotal)
    # dihedral_mutator.stepsize *= 0.95
    # Print.as_xyz(st, ostream = file_xyz, title = "Step: $step")
end
my_callback = Common.CallbackObject(1, callback)

#9. DEFINE THE DRIVER PARAMETERS AND RUN THE SIMULATION
mc_driver = Drivers.MonteCarlo.MonteCarloDriver(my_sampler, my_evaluator!, 1.0, 10)
Drivers.MonteCarlo.run!(state, mc_driver, my_callback)

close(file_xyz)
exit(1)
# ------------------------------------------------------------------------------------




# ---
#Add dummy dihedral for first C_alpha (To correctly find the first C alpha)
insert!(dihedrals, 1, Mutators.Dihedral.NewDihedral(-1, 2, 3, length(residues[1].atoms) - 1, residues[1].atoms[4:end], residues[1], "PHI"))
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