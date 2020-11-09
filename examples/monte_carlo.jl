# Load resources
include("../src/ProtoSyn.jl")

using .ProtoSyn
using .ProtoSyn.Calculators
using .ProtoSyn.Peptides
using .ProtoSyn.Builder
using .ProtoSyn.Units
println(" | ProtoSyn loaded successfully.")

using Printf
using Random

let
    seed = convert(Int, floor(rand() * 100000))
    # seed = 46747
    Random.seed!(seed)

    # Start algorithm
    res_lib = grammar();
    pose = Peptides.build(res_lib, seq"MGSWAEFKQRLAAIKTRLQALGGSEAELAAFEKEIAAFESELQAY");
    Peptides.setss!(pose, SecondaryStructure[:helix], rid"1:20");
    Peptides.setss!(pose, SecondaryStructure[:helix], rid"27:44");
    sync!(pose)

    energy_function = Calculators.get_energy_function()

    movable_residues = 20:27
    n_starts         = 10
    n_cycles         = 5
    n_steps          = 10000
    init_temp        = 0.05
    print_every      = 100
    init_step_size   = 1.0
    accepted         = 0
    e                = energy_function(pose)
    init_e           = e
    be               = e
    results          = Dict(1 => e)
    saved_state      = copy(pose.state)
    init_state       = copy(pose.state)
    best_state       = copy(pose.state)
    io               = open("../teste.pdb", "w"); ProtoSyn.write(io, pose); close(io);
    io               = open("../teste_run.dat", "w"); close(io)


    # Loop algorithm
    for start in 1:n_starts
        for cycle in 1:n_cycles
            # start_temp = - (((cycle - 1)^2 * init_temp) / (n_cycles^2)) + init_temp # quadratic
            start_temp = - (((cycle - 1) * init_temp) / n_cycles) + init_temp # linear
            # start_step_size = - (((cycle - 1)^2 * init_step_size) / (n_cycles^2)) + init_step_size # quadratic
            start_step_size = - (((cycle - 1) * init_step_size) / n_cycles) + init_step_size # linear
            
            for step in 1:n_steps    

                # temperature = - (((step - 1)^2 * start_temp) / (n_steps^2)) + start_temp # quadratic
                temperature = - (((step - 1) * start_temp) / n_steps) + start_temp
                # step_size = - (((step - 1)^2 * start_step_size) / (n_steps^2)) + start_step_size # quadratic
                step_size = - (((step - 1) * start_step_size) / n_steps) + start_step_size # linear
            
                r = rand(pose.graph.items[1].items[movable_residues])
                phi_psi = rand([Peptides.Dihedral.phi, Peptides.Dihedral.psi])
                val = randn() * π * step_size # in radians
                Peptides.rotate_dihedral!(pose.state, r, phi_psi, val)
                ne = energy_function(pose)
                n = rand()

                if (ne < e) || (n < exp((e - ne)/ temperature))
                    e = ne
                    saved_state = copy(pose.state)
                    accepted += 1
                else
                    pose.state = copy(saved_state)
                end

                if ne < be
                    best_state = copy(pose.state)
                    be = ne
                end

                if step % print_every == 0
                    status = @sprintf("Start: %2d | Cycle: %2d | Step: %4d | Energy: %12.5f | Acceptance Rate: %6.4f (%-4d) | Temperature: %6.4f | Step size: %6.4f\n", start, cycle, step, e, accepted/step, accepted, temperature, step_size) 
                    print(status)
                    # io = open("../teste.pdb", "a"); ProtoSyn.write(io, pose); close(io);
                    io = open("../teste_run.dat", "a")
                    write(io, status)
                    close(io);
                end
            end
            results[start + start * (cycle + cycle * ceil(n_steps / print_every))] = be
            println("\n | Results: $results")

            # Recover best state
            pose.state  = copy(best_state)
            saved_state = copy(best_state)
            accepted    = 0
            e  = be
            io = open("../teste.pdb", "a"); ProtoSyn.write(io, pose); close(io);
        end

    pose.state  = copy(init_state)
    saved_state = copy(init_state)
    best_state  = copy(init_state)
    accepted    = 0
    e           = init_e
    be          = init_e
    end
    println(" | Seed: $seed")
end