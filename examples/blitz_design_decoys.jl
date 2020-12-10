# In this example, we will attempt to launch N number of decoys of a simulation,
# evaluating the best decoy at the end. We can do this in serial mode (1
# simulation at a time), or in parallel, where multiple processors (even in
# multiple machines) do 1 simulation each one (at the same time) and the results
# are compared at the end, therefore speeding up the calculation time. In a
# single machine (running Ubuntu), the number of available threads can be found
# with the command: cat /proc/cpuinfo | grep processor | wc -l

using Distributed

# We must first define the number of workers we want to make available to
# ProtoSyn, usually using the `addprocs` function.
addprocs(4)

# Loaded packages are loaded in all worker processes
using ProtoSyn
using ProtoSyn.Peptides
using Printf

# The function that actually does the work needs to be loaded in all worker
# processes, by using the @everywhere macro. This function is just a condensed
# version of the blitz_design.jl example.
@everywhere function work(jobs, results)
    while true
        job_id = take!(jobs)
        # T = Float64
        # res_lib = grammar(T)
        # rot_lib = Peptides.Rotamers.load_dunbrack(T)
        # pose = Peptides.load("examples/data/mdC.pdb")
        # energy_function = ProtoSyn.Common.get_default_energy_function()
        # design_selection = (5.0:rn"CBZ") & !(rn"CBZ" | rn"PRO")
        # design_selection = ProtoSyn.promote(design_selection, Residue)
        # p                = 1/count(design_selection(pose))
        # design_mutator   = Peptides.Mutators.DesignMutator(p, res_lib, design_selection)
        # filter!(aa -> aa != 'P', design_mutator.searchable_aminoacids)
        # rotamer_blitz = ProtoSyn.Drivers.RotamerBlitz(energy_function, rot_lib, 3, 1)
        # function my_custom_sampler!(pose::Pose)
        #     pre_design_sequence  = Peptides.get_sequence(pose)
        #     design_mutator(pose)
        #     post_design_sequence = Peptides.get_sequence(pose)
        #     blitz_selection = !TrueSelection{Residue}()
        #     for i in 1:length(pre_design_sequence)
        #         if pre_design_sequence[i] != post_design_sequence[i]
        #             blitz_selection = blitz_selection | SerialSelection{Residue}(i, :id)
        #         end
        #     end
        #     rotamer_blitz.selection = blitz_selection
        #     rotamer_blitz(pose)
        # end
        # monte_carlo = ProtoSyn.Drivers.MonteCarlo(
        #     energy_function,
        #     my_custom_sampler!,
        #     nothing,
        #     100,
        #     ProtoSyn.Drivers.get_linear_quench(0.1, 100)
        # )
        # monte_carlo(pose)
        # put!(results, pose)
    end
end

# A worker will take a job, finish it, and then look for another job, until no
# more jobs are available. The list of jobs is just a way for us to share work
# equitably between workers (a worker might be able to perform 2 simulations
# while another worker does a single one, because the simulations were smaller
# or the resources available were different). When no more jobs exist on the
# `jobs` channel, the workers stop. N is the number of simulations requested.
N    = 8
jobs = RemoteChannel(() -> Channel{Int32}(N))

# The results are `put!` in a RemoteChannel. This is just a shared version of
# Channel, which in turn acts as a "pipe" that receives and outputs a certain
# type of variable (Pose, is this case). Each worker, one its simulation is
# finished, will `put!` the resulting pose in this channel, for the master
# process to collect and analyze.
results = RemoteChannel(() -> Channel{Pose}(N))

# We can now give work to the workers. `make_jobs` simply creates a list of
# Int32 instances (like job cards).
function make_jobs(n)
    for i in 1:n
        put!(jobs, i)
    end
end

make_jobs(N)

# `remote_do` automatically spreads the tasks by the available workers.
for p in workers()
    remote_do(work, p, jobs, results)
end