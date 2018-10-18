
using ProtoSyn
using Printf

dihd_json = Aux.read_JSON("data/1ctf_mc_top.json")
dihedrals, residues = Mutators.Dihedral.load_topology(dihd_json);

state = Common.load_from_pdb("data/mol.pdb")
state.energy = Common.Energy();

dihedral_mutator = Mutators.Dihedral.DihedralMutator(
    dihedrals   # list of dihedrals
    , 0.1       # single dihedral mutation probability
    , randn     # angle sampler
    , 1.0       # stepsize
)
my_sampler(st::Common.State) = Mutators.Dihedral.run!(st, dihedral_mutator);

function my_evaluator(st::Common.State, do_forces::Bool)
    st.energy.eTotal
end;

mc_driver = Drivers.MonteCarlo.MonteCarloDriver(
    my_sampler      # sampler
    , my_evaluator  # evaluator
    , 1.0           # temperature
    , 1000          # nsteps
);

function callback(step::Int64, st::Common.State, dr::Drivers.MonteCarlo.MonteCarloDriver, ac_ratio::Float64)
    write(stdout, @sprintf "(MC) Step: %4d | Energy: %9.4f\n" step state.energy.eTotal)
    dihedral_mutator.stepsize *= 0.95
end

my_callback = Common.CallbackObject(
    100          # calling frequency
    , callback  # the actual callback function
);



cb = @callback 10 function(step::Int64, st::Common.State, dr::Drivers.MonteCarlo.MonteCarloDriver, naccept::Int64)
    # do something
    println("do something")
end

Drivers.MonteCarlo.run!(state, mc_driver, my_callback)
