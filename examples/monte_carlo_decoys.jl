using Distributed

n_workers = 4
addprocs(n_workers)

@everywhere using ProtoSyn
@everywhere using ProtoSyn
@everywhere using ProtoSyn.Peptides
@everywhere using ProtoSyn.Calculators
@everywhere using Printf

@everywhere begin
    
    T        = Float64
    res_lib  = Peptides.grammar(T)
    sequence = seq"MGSWAEFKQRLAAIKTRLQALGGSEAELAAFEKEIAAFESELQAYKGKGNPEVEALRKEAAAIRDELQAYRHN"

    pose = Peptides.build(res_lib, sequence);

    Peptides.setss!(pose, SecondaryStructure[:helix], rid"1:20");
    Peptides.setss!(pose, SecondaryStructure[:helix], rid"27:44");
    Peptides.setss!(pose, SecondaryStructure[:helix], rid"50:end");

    Peptides.remove_sidechains!(pose)
    sync!(pose)

    energy_function = ProtoSyn.Common.default_energy_function()
    GC_frequency = round(Int16, 6)
    energy_function = Calculators.EnergyFunction(
        Dict{Calculators.EnergyFunctionComponent, T}(
            Calculators.TorchANI.torchani_model => T(1.0),
            Peptides.Calculators.Caterpillar.solvation_energy => T(0.05),
            Calculators.Restraints.bond_distance_restraint => T(0.05)),
        GC_frequency)

    selection          = (rid"21:26" | rid"45:49") & an"C$|N$"r
    p_mut              = 1/(count(selection(pose)))
    dihedral_mutator   = ProtoSyn.Mutators.DihedralMutator(
        randn, p_mut, 0.5, selection)

    selection          = (rid"21:26" | rid"45:49") & an"CA" & !rn"PRO"
    n                  = count(selection(pose))
    p_mut              = 0.5/(n*(n-1))
    crankshaft_mutator = ProtoSyn.Mutators.CrankshaftMutator(
        randn, p_mut, 0.05, selection, !(an"^CA$|^N$|^C$|^H$|^O$"r))

    compound_driver = ProtoSyn.Drivers.CompoundDriver(
        [crankshaft_mutator, dihedral_mutator])

    # function callback_function(pose::Pose, driver_state::ProtoSyn.Drivers.DriverState, other...)
    #     @printf("STEP %-6d | E= %-10.4f | AR= %-5.1f%%  | T= %-7.4f | MOVS: C-%d D-%d\n", driver_state.step, pose.state.e[:Total], (driver_state.acceptance_count/driver_state.step)*100, driver_state.temperature, other[1][1], other[1][2])
    #     driver_state.step == 0 && return
    #     ProtoSyn.append(pose, "monte_carlo.pdb", model = driver_state.step)
    # end

    # callback    = ProtoSyn.Drivers.Callback(callback_function, 10)
    n_steps     = 100_000
    monte_carlo = ProtoSyn.Drivers.MonteCarlo(
        energy_function,
        compound_driver,
        nothing,
        n_steps,
        ProtoSyn.Drivers.get_linear_quench(0.4, n_steps))

    function start_simulation(job_cards::RemoteChannel, results::RemoteChannel)
        while true
            job = take!(job_cards)
            pose = Peptides.build(res_lib, sequence);

            Peptides.setss!(pose, SecondaryStructure[:helix], rid"1:20");
            Peptides.setss!(pose, SecondaryStructure[:helix], rid"27:44");
            Peptides.setss!(pose, SecondaryStructure[:helix], rid"50:end");
        
            Peptides.remove_sidechains!(pose)
            sync!(pose)

            monte_carlo(pose)
            put!(results, pose)
        end
    end
end

N = 20

job_queue = RemoteChannel(() -> Channel{Int16}(N))
for job_n in 1:N
    put!(job_queue, Int16(job_n))
end

results_queue = RemoteChannel(() -> Channel{Any}(N))

for p in workers()
    remote_do(start_simulation, p, job_queue, results_queue)
end

begin
    local n = N
    results  = []
    while n > 0
        result = take!(results_queue)
        push!(results, result)
        n = n - 1
    end
end

ProtoSyn.write(results[1], "blitz_design_decoys.pdb")
for index in 2:length(results)
    ProtoSyn.append(results[index], "blitz_design_decoys.pdb", model = index)
end