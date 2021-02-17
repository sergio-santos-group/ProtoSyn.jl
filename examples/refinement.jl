using ProtoSyn
using ProtoSyn.Peptides
using ProtoSyn.Builder
using Printf

res_lib  = ProtoSyn.Peptides.grammar()

sequence = seq"MGSWAEFKQRLAAIKTRLQALGA"
pose     = ProtoSyn.Peptides.build(res_lib, sequence)

begin
    pose = ProtoSyn.Peptides.load("monte_carlo_1.pdb")
    ProtoSyn.Peptides.add_sidechains!(pose, res_lib)
end
# ProtoSyn.write(pose, "refinement.pdb")
energy_function = ProtoSyn.Common.default_energy_function()
cb = ProtoSyn.Common.default_energy_step_frame_callback(10, "refinement.jl")
sd = ProtoSyn.Drivers.SteepestDescent(energy_function, cb, 1000, 0.001, 0.1)
sd(pose)



α_sol   = 0.01
α_cmap  = 0.001
T_i     = 0.1
n_steps = 100

energy_function = ProtoSyn.Common.default_energy_function()
cm_restraints   = Peptides.Calculators.Restraints.load_contact_map("contact_map_example.txt")
energy_function.components[ProtoSyn.Peptides.Calculators.Caterpillar.solvation_energy] = α_sol
energy_function.components[cm_restraints] = α_cmap
energy_function(pose)

sd_energy_function = copy(energy_function)
sd_energy_function.components[cm_restraints] = 0.0
cb = ProtoSyn.Common.default_energy_step_callback(500)
sd = ProtoSyn.Drivers.SteepestDescent(sd_energy_function, cb, 1000, 0.001, 0.1)

brm = ProtoSyn.Mutators.BlockRotMutator(ProtoSyn.rand_vector_in_sphere, randn, ProtoSyn.center_of_mass, 0.3, 0.3, [rid"1:20", rid"27:44", rid"50:end"], nothing)
btm = ProtoSyn.Mutators.BlockTransMutator(ProtoSyn.rand_vector_in_sphere, 0.3, 0.3, [rid"1:20", rid"27:44", rid"50:end"], sd)

compound_driver = ProtoSyn.Drivers.CompoundDriver([brm, btm])

function callback_function(pose::Pose, driver_state::ProtoSyn.Drivers.DriverState)
    s = @sprintf(" INNER STEP %-6d | E(total)= %-10.4f | AR= %-5.1f%%  | T= %-7.4f | E(ani)= %-10.4f | E(sol)= %-10.4f | E(cmp)= %-10.4f | E(cls)= %-10.4f\n", driver_state.step, pose.state.e[:Total], (driver_state.acceptance_count/driver_state.step)*100, driver_state.temperature, pose.state.e[:TorchANI_ML_Model], pose.state.e[:Caterpillar_Solvation], pose.state.e[:Contact_Map_Restraint], pose.state.e[:Clash_Restraint])
    print(s)
    open("monte_carlo.dat", "a") do io
        write(io, s)
    end
    driver_state.step == 0 && return
    ProtoSyn.append(pose, "monte_carlo.pdb")
end

callback    = ProtoSyn.Drivers.Callback(callback_function, 1)
monte_carlo = ProtoSyn.Drivers.MonteCarlo(
    energy_function,
    compound_driver,
    callback,
    n_steps,
    ProtoSyn.Drivers.get_linear_quench(T_i, n_steps))

ProtoSyn.write(pose, "monte_carlo.pdb")
monte_carlo(pose)