# Load resources
include("../src/ProtoSyn.jl")

using .ProtoSyn
using .ProtoSyn.Peptides
using .ProtoSyn.Builder
using .ProtoSyn.Units
println(" | ProtoSyn loaded successfully.")

using Printf
using Random
using PyCall

let
    seed = convert(Int, floor(rand() * 100000))
    # seed = 46747
    Random.seed!(seed)

    torch    = pyimport("torch")
    torchani = pyimport("torchani")

    if torch.cuda.is_available()
        device = torch.device("cuda")
    else
        device = torch.device("cpu")
    end

    model = torchani.models.ANI2x(periodic_table_index = true).to(device)
    println(" | TorchANI model loaded successfully.")

    # Define functions
    function calc_energy(pose::Pose, model)
        c = get_cartesian_matrix(pose)
        coordinates = torch.tensor([c], requires_grad = true, device = device).float()
        ani_energy = float(model((species, coordinates))[2])

        return ani_energy
    end

    function model2(t::Tuple)
    
        m1 = model.species_converter(t)
        m2 = model.aev_computer(m1)
        m3 = get(model.neural_networks, 0)(m2)
        return m3
    end


    # Start algorithm
    res_lib = grammar();
    pose = Peptides.build(res_lib, seq"MGSWAEFKQRLAAIKTRLQALGGSEAELAAFEKEIAAFESELQAY");
    Peptides.setss!(pose, SecondaryStructure[:helix], rid"1:20");
    Peptides.setss!(pose, SecondaryStructure[:helix], rid"27:44");

    _model           = model
    movable_residues = 20:27
    s                = get_ani_species(pose)
    species          = torch.tensor([s], device = device)
    n_starts         = 10
    n_cycles         = 5
    n_steps          = 3000
    init_temp        = 0.1
    print_every      = 100
    init_step_size   = 1.0
    accepted         = 0
    e                = calc_energy(pose, _model)
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
                ne = calc_energy(pose, _model)
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
                    io = open("../teste.pdb", "a"); ProtoSyn.write(io, pose); close(io);
                    io = open("../teste_run.dat", "a")
                    write(io, status)
                    close(io);
                end
            end
            results[start + start * (cycle + cycle * ceil(n_steps / print_every))] = be
            # println("\n | Results: $results")

            pose.state  = copy(best_state)
            saved_state = copy(best_state)
            accepted    = 0
            # e           = calc_energy(pose, _model)
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