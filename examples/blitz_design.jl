# In this example we will try to create a monte carlo simulation for design
# around a given ligand. In each step of the simulation we want to change N
# residues and then do a blitz rotamer search on the residues surrouding the
# chance, stabilizing it.

using Dates

start_time = now()

using ProtoSyn
using ProtoSyn.Peptides
using Printf

# 1) Load the default residue library
T = Float64
res_lib = grammar(T)

# 2) Load the default Dunbrack rotamer library
rot_lib = Peptides.Rotamers.load_dunbrack(T)

# 3) Load the initial pose
pose = Peptides.load("data/mdC.pdb")

# 5) Load the default energy function
energy_function = ProtoSyn.Common.get_default_energy_function()

# 4) Define the DesignMutator (Here we will be ignoring the prolines, and will
# target the residues within 5 Å of the ligand "CBZ").
design_selection = (5.0:rn"CBZ") & !(rn"CBZ" | rn"PRO")
design_selection = ProtoSyn.promote(design_selection, Residue)
p                = 1/count(design_selection(pose))
design_mutator   = Peptides.Mutators.DesignMutator(p, res_lib, design_selection)
# Remove proline from searchable aminoacids
filter!(aa -> aa != 'P', design_mutator.searchable_aminoacids)

# 5) Define the RotamerBlitz driver (Here we will consider the 3 most likely
# rotamers, with just 1 pass).
rotamer_blitz = ProtoSyn.Drivers.RotamerBlitz(energy_function, rot_lib, 3, 1)

# 6) Define the custom sampler
# For this example, we want to identify the changed residues (each step) and
# perform a blitz rotamer search localized exclusevly on the surrouding
# residues. Therefore, this extra step of identifying mutated residues needs to
# be explicitly employed, and so we can't use a CompoundMutator, but instead
# need to define a custom function. This function should receive and change
# "in-place" a pose. In order to tell the RotamerBlitz driver what residue to
# consider, we must update it's `selection` parameter.
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

# 7) Create the Monte Carlo driver (Here we will do a linear temperature
# quenching over 100 simulation steps, with a callback every 5 saving the
# structure to a file).
function callback(pose::Pose, driver_state::ProtoSyn.Drivers.DriverState)
    @printf("STEP %-5d E= %-10.4f AR= %-5.2f%% T= -%7.5f\n", driver_state.step, pose.state.e[:Total], (driver_state.acceptance_count/driver_state.step)*100, driver_state.temperature)
    ProtoSyn.append(pose, "../blitz_design.pdb", model = driver_state.step)
end

mc_callback = ProtoSyn.Drivers.Callback(callback, 5)
n_steps = 1000
monte_carlo = ProtoSyn.Drivers.MonteCarlo(
    energy_function,
    my_custom_sampler!,
    nothing,
    n_steps,
    ProtoSyn.Drivers.get_linear_quench(0.1, n_steps)
)

# ProtoSyn.write(pose, "../blitz_design.pdb")

# 8) Run the simulation
N = 20
for i in 1:N
    pose = Peptides.load("data/mdC.pdb")
    monte_carlo(pose)
    if i == 1
        ProtoSyn.write(pose, "blitz_design.pdb")
    else
        ProtoSyn.append(pose, "blitz_design.pdb", model = i)
    end
end

elapsed_time = Dates.canonicalize(Dates.CompoundPeriod(now() - start_time))
println("| All tasks completed in $elapsed_time.\n")