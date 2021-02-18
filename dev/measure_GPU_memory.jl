using ProtoSyn
using Printf

energy_function = ProtoSyn.Calculators.EnergyFunction(
    Dict(ProtoSyn.Calculators.TorchANI.torchani_model => Float64(1.0)))

D        = 1.5
R        = 10
filename = "gpu_results.csv"
N_steps  = 1000

pose = ProtoSyn.Materials.Lattices.primitive()
ProtoSyn.symexp!(pose, [R, R, R], [D, D, D])

energy_function(pose, update_forces = true)
# ProtoSyn.Calculators.TorchANI.calc_torchani_model_xmlrpc(pose, update_forces = true)

open(filename, "a") do file_out
    # write(file_out, @sprintf("%10s %-10s %-12s %-10s %-20s\n", "N atoms", "GC", "GC (Forces)", "XML-RPC", "XML-RPC (Forces)"))

    println("GC - No forces")
    start = time()
    for step in 1:N_steps
        energy_function(pose, update_forces = false)
    end
    gc_nf = time() - start
    GC.gc()
    println(" -> Allocation: $(ProtoSyn.gpu_allocation())")

    println("GC - forces")
    start = time()
    for step in 1:N_steps
        energy_function(pose, update_forces = true)
    end
    gc_f = time() - start
    GC.gc()
    println(" -> Allocation: $(ProtoSyn.gpu_allocation())")

    # println("XMLRPC - No forces")
    # start = time()
    # for step in 1:N_steps
    #     ProtoSyn.Calculators.TorchANI.calc_torchani_model_xmlrpc(pose, update_forces = false)
    # end
    xml_nf = 0.0
    # GC.gc()
    # println(" -> Allocation: $(ProtoSyn.gpu_allocation())")

    # println("XMLRPC - Forces")
    # start = time()
    # for step in 1:N_steps
    #     ProtoSyn.Calculators.TorchANI.calc_torchani_model_xmlrpc(pose, update_forces = true)
    # end
    xml_f = 0.0
    # GC.gc()
    # println(" -> Allocation: $(ProtoSyn.gpu_allocation())")

    s = @sprintf("%10d %-10.3f %-12.3f %-10.3f %-20.3f\n", pose.state.size, gc_nf, gc_f, xml_nf, xml_f)
    write(file_out, s)
end

ProtoSyn.Calculators.TorchANI.stop_torchANI_server()