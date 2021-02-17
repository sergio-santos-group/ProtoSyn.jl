using ProtoSyn
using Printf
using CUDA

energy_function = ProtoSyn.Common.default_energy_function()
energy_function.components[ProtoSyn.Calculators.Restraints.clash_restraint] = 0.0
energy_function.components[ProtoSyn.Calculators.Restraints.bond_distance_restraint] = 0.0
energy_function.components[ProtoSyn.Peptides.Calculators.Caterpillar.solvation_energy] = 0.0

D = 1.5
R = 4
pose = ProtoSyn.Materials.Lattices.primitive()
ProtoSyn.symexp!(pose, [R, R, R], [D, D, D])
energy_function(pose)

filename = "gpu_results.csv"
N_atoms  = pose.state.size
N_steps  = 10000

function gpu_allocation()
    return (CUDA.total_memory() - CUDA.available_memory()) / CUDA.total_memory()
end

open(filename, "w") do file_out
    write(file_out, @sprintf("N atoms: %10d\n", N_atoms))
    write(file_out, @sprintf("%10s %12s %10s %12s\n", "Step #", "Allocation %", "Elapsed s", "Energy e.u."))
    write(file_out, @sprintf("%10d %12.3f %10.3f %12s\n", 0, gpu_allocation(), 0.0, "None"))

    start = time()
    for step in 1:N_steps
        e = energy_function(pose)
        s = @sprintf("%10d %12.3f %10.3f %12.3f\n", step, gpu_allocation(), time() - start, e)
        write(file_out, s)
        println(s)
    end
end