# In this example, we will attempt to launch N number of decoys of a simulation,
# evaluating the best decoy at the end. We can do this in serial mode (1
# simulation at a time), or in parallel, where multiple processors (even in
# multiple machines) do 1 simulation each one (at the same time) and the results
# are compared at the end, therefore speeding up the calculation time. In a
# single machine (running Ubuntu), the number of available threads can be found
# with the command: cat /proc/cpuinfo | grep processor | wc -l
# More information regarding parallel computing in Julia environments can be
# found at https://docs.julialang.org/en/v1/manual/distributed-computing/

using Dates

start_time = now()

using Distributed

# Define the number of workers we want to make available to ProtoSyn
n_workers = 4
addprocs(n_workers)

# It's necessary to load the ProtoSyn package in all workers. Therefore, this
# can only be done after `addprocs`. The Distributed package loaded earlier is a
# special caso which is automatically loaded in all workers.
@everywhere using ProtoSyn
@everywhere using ProtoSyn.Peptides
@everywhere using ProtoSyn.Builder
@everywhere using ProtoSyn.Calculators
@everywhere using Printf

# We will divide our simulation in two parts: the 'setup' and the 'simulation'.
# The setup stage should run only once in each worker, setting up the necessary
# variables, while the simulation part should run N times, where N is the number
# of replicas desired.
# (A) Setup stage: The variables loaded in this stage are further explained in
# the blitz_design.jl example.
@everywhere begin
    T = Float64
    res_lib = grammar(T)
    rot_lib = Peptides.Rotamers.load_dunbrack(T)
    pose    = Peptides.load("data/mdC.pdb") # Only to define p_mut

    # For a distributed parallel computing application, we must consider the
    # avaliability of resources. In this case, TorchANI model uses the GPU,
    # caching a certain amount of memory for the model. Each worker will need
    # a part of this memory, and so the number of local workers might be limited
    # due to lack of GPU memory (depending on the size of the simulated system).
    # In order to help reduce the memory fingerprint, we can increase the
    # garbage collection frequency. However, this process can be costly (in
    # terms of time), and so a controlled equilibrium in this parameter should
    # always be a concern in parallel computing applications.
    GC_frequency = round(Int16, 50/8)
    energy_function = Calculators.EnergyFunction(
        Dict{Calculators.EnergyFunctionComponent, T}(
            Calculators.TorchANI.torchani_model => T(1.0)),
        GC_frequency)

    design_selection = (5.0:rn"CBZ") & !(rn"CBZ" | rn"PRO")
    design_selection = ProtoSyn.promote(design_selection, Residue)
    p                = 1/count(design_selection(pose))
    design_mutator   = Peptides.Mutators.DesignMutator(p, res_lib, design_selection)
    filter!(aa -> aa != 'P', design_mutator.searchable_aminoacids)
    rotamer_blitz = ProtoSyn.Drivers.RotamerBlitz(energy_function, rot_lib, 3, 1)
    
    function my_custom_sampler!(pose::Pose)
        pre_design_sequence  = Peptides.get_sequence(pose)
        design_mutator(pose)
        post_design_sequence = Peptides.get_sequence(pose)
        
        blitz_selection = !TrueSelection{Residue}()
        for i in 1:length(pre_design_sequence)
            if pre_design_sequence[i] != post_design_sequence[i]
                blitz_selection = blitz_selection | SerialSelection{Residue}(i, :id)
            end
        end
    
        rotamer_blitz.selection = blitz_selection
        rotamer_blitz(pose)
    end

    function callback(pose::Pose, driver_state::ProtoSyn.Drivers.DriverState)
        @printf("STEP %-5d E= %-10.4f AR= %-5.2f%% T= -%7.5f\n", driver_state.step, pose.state.e[:Total], (driver_state.acceptance_count/driver_state.step)*100, driver_state.temperature)
    end
    
    mc_callback = ProtoSyn.Drivers.Callback(callback, 5)

    n_steps = 1000
    monte_carlo = ProtoSyn.Drivers.MonteCarlo(
        energy_function,
        my_custom_sampler!,
        mc_callback,
        n_steps,
        ProtoSyn.Drivers.get_linear_quench(0.1, n_steps)
    )

    # Consume the job queue until exhausted
    function start_simulation(job_cards::RemoteChannel, results::RemoteChannel)
        while true
            job = take!(job_cards)
            pose = Peptides.load("data/mdC.pdb")
            monte_carlo(pose)
            put!(results, pose)
        end
    end
end

# (B) Simulation Stage: There are different approaches to spawn N jobs in the
# workers. For this example, we will employ a queue consumption method - all the
# requested jobs are placed in a queue. Once free, workers will consume a job
# from the queue, until no more jobs are available. Since the workers can be in
# a remote machine, a RemoteChannel is employed, rather than a simple Channel.
# Set the number of replicas
N = 20

# Create a job queue with N job cards. A job card is simply an Int16, could be
# anything. The objetive is to take job_cards from the queue until the requested
# number of replicas is performed. This allows for equitable sharing of work
# between workers (1 simulation might take longer than the other, the rest of
# workers can take new jobs without synching among themselves).
job_queue = RemoteChannel(() -> Channel{Int16}(N))
for job_n in 1:N
    put!(job_queue, Int16(job_n))
end

# # Create results_queue where replica results will be stored
results_queue = RemoteChannel(() -> Channel{Any}(N))

# Spawn the workers. The `start_simulation`function has a `while true` loop
# which will run until no more jobs exist in the queue.
for p in workers()
    remote_do(start_simulation, p, job_queue, results_queue)
end

# After starting the job_queue consumption, wait for all replicas to be
# completed, synching the workers
begin
    local n = N
    results  = []
    while n > 0
        result = take!(results_queue)
        push!(results, result)
        n = n - 1
    end
end

# Write results to multiple model PDB file
ProtoSyn.write(results[1], "blitz_design_decoys.pdb")
for index in 2:length(results)
    ProtoSyn.append(results[index], "blitz_design_decoys.pdb", model = index)
end

elapsed_time = Dates.canonicalize(Dates.CompoundPeriod(now() - start_time))
println("| All tasks completed in $elapsed_time.\n")

# Note: In a benchmark test, the same simulation with 1000 steps and 20 replicas
# took 26 minutes in sequential mode, but only 11 minutes in parallel.