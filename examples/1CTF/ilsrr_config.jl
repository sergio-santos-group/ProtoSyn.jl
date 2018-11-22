# ------------------
# SET CONFIGURATION:
# ------------------

# 1. Define the algorithm length:
const n_outer_steps            = 1000   # Number of inner cycles to perform;
const n_inner_steps            = 250    # Number of steps on each inner cycle;
const init_min_n_steps         = 100000  # Number of steps for Steepest Descent minimization (initial long minimization);
const min_n_steps              = 10     # Number of steps for Steepest Descent minimization (every step of inner cycle);

# 2. Define output frequency:
const print_sts_every_ic       = 10     # Print system status every n steps;
const print_str_every_ic       = 10     # Print structure to files every n steps;
const print_sts_every_min      = 250     # Print status every n steps in Steepest Descent minimization;

# 3. Define the step size:
const min_step_size            = 0.01   # Minimum step_size allowed;
const updt_step_size_every     = 0      # Update step_size every n steps;
dihedral_step_size             = π/2    # (Initial) Step size for dihedral movements;
crankshaft_step_size           = π/2    # (Initial) Step size for crankshaft movements;
const dihedral_step_size_prt   = π/2    # Step size for dihedral movements in the perturbator;
const crankshaft_step_size_prt = π/2    # Step size for crankshaft movements in the perturbator;

# 4. Define the mutation probabilities:
const dihedral_p_mut           = 0.025  # Mutation probability for dihedral movements;
const dihedral_p_mut_prt       = 0.25   # Mutation probability for dihedral movements in the perturbator;
const crankshaft_p_mut         = 0.007  # Mutation probability for crankshaft movements;
const crankshaft_p_mut_prt     = 0.07   # Mutation probability for crankshaft movements in the perturbator;

# 5. Define the system (Initial) temperature:
const inner_temperature        = 300.0  # Temperature for Metropolis acceptance in inner cycles (0.0 performs a Random Search where only better structures are accepted);
const outer_temperature        = 1000.0 # Temperature for Metropolis acceptance in outer cycles. The bigger, the easier to accept worse configurations as new homebases;

# 6. Define the target acceptance ratio and buffer zone:
const acceptance_ratio         = 0.2    # Target acceptance ration of new structures. Step_size will be updated to try and achieve this value;
const ar_buffer_zone           = 0.01   # Buffer zone limit around the target acceptance_ratio. Step_size will not be updated if it's currently inside the buffer zone;

# 7. Define the molecule secondary structure
const ss                       = "CEEEEEEECCCCHHHHHHHHHHHHCCCHHHHHHHHHCCCEEEEEEECHHHHHHHHHHHHHHCCEEEEC"

# 8. Define the input files:
const input_pdb                = "data/1ctf_good_rmsd.pdb"
const input_contact_map        = "data/contact_map_gremlin.txt"
const input_amber_json         = "data/1ctf_amber_top.json"

# 9. Define the output destination:
const log_destination          = stdout
const xyz_destination          = open("out/trajectory.pdb", "w")
const best_destination         = open("out/best_trajectory.pdb", "w")
const out_destination          = open("out/data.out", "w")
const outer_best_destination   = open("out/outer_best_destination.pdb", "w")
const status_destination       = [log_destination, out_destination]