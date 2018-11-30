# ------------------
# SET CONFIGURATION:
# ------------------

# 1. Define the algorithm length:
const n_search_refine_cycles   = 10     # Number of cycles between Random Search / Refinement stages;
const n_i_search_steps         = 10000  # Number of steps on each Monte Carlo Search inner cycle;
const n_o_search_steps         = 100    # Number of ILSRR outer cycles to perform;
const n_init_min_steps         = 10000  # Number of steps for Steepest Descent minimization (initial long minimization);
const n_refine_min_steps       = 5000   # Number of steps for Steepest Descent minimization (every step of Monte Carlo Refinement);
const n_loop_closer_steps      = 20000  # Number of steps for Steepest Descent in loop closure


# 2. Define output frequency:
const print_sts_every_ic       = 1      # Print system status every n steps (Inner cycle);
const print_sts_every_oc       = 1      # Print system status every n steps (Outer cycle);
const print_sts_every_min      = 1000   # Print status every n steps in Steepest Descent minimization;
const print_str_every_ic       = 50     # Print structure to files every n steps;
const print_str_every_min      = 0      # Print structure to files every n steps in Steepest Descent minimization;


# 3. Define the step size:
const min_step_size            = 0.01   # Minimum step_size allowed;
const soft_dh_step_size        = π/4    # Step size for dihedral movements in the Monte Carlo Refinement;
reg_dh_step_size               = π/2    # (Initial) Step size for dihedral movements in the Monte Carlo Search (ILSRR inner cycle);
const hard_dh_step_size        = π/2    # Step size for dihedral movements in the ILSRR perturbator;

const soft_cs_step_szie        = π/4    # Step size for crankshaft movements in the Monte Carlo Refinement;
reg_cs_step_size               = π/2    # (Initial) Step size for crankshaft movements in the Monte Carlo Search (ILSRR inner cycle);
const hard_cs_step_size        = π/2    # Step size for crankshaft movements in the ILSRR perturbator;

const soft_br_step_size        = π/8    # Step size for blockrot movements in the Monte Carlo Refinement;


# 4. Define convergence criteria for Steepest Descent:
const max_step                 = 1e-4
const f_tol                    = 50.0


# 5. Define the mutation probabilities:
const soft_dh_p_mut            = 0.025  # Mutation probability for dihedral movements in the Monte Carlo Refinement;
const reg_dh_p_mut             = 0.025  # Mutation probability for dihedral movements in the Monte Carlo Search (ILSRR inner cycle);
const hard_dh_p_mut            = 0.25   # Mutation probability for dihedral movements in the ILSRR perturbator;

const soft_cs_p_mut            = 0.007  # Mutation probability for crankshaft movements in the Monte Carlo Refinement;
const reg_cs_p_mut             = 0.007  # Mutation probability for crankshaft movements in the Monte Carlo Search (ILSRR inner cycle);
const hard_cs_p_mut            = 0.07   # Mutation probability for crankshaft movements in the ILSRR perturbator;

const soft_br_p_mut            = 0.15   # Mutation probability for blockrot movements in the Monte Carlo Refinement;


# 6. Define the system temperature(s):
i_search_temp_init             = 300.0  # Temperature for Metropolis acceptance in inner cycles (0.0 performs a Random Search where only better structures are accepted);
const o_search_temp_static     = 300.0  # Temperature for ILSRR acceptance in outer cycles (0.0 performs a Random Search where only better structures are accepted);
const refine_temp_static       = 300.0  # Temperature for Metropolois acceptance in refinement (0.0 performs a Random Search where only better structures are accepted);


# 7. Define the target acceptance ratio for all Monte Carlo Drivers and respective buffer zone:
const acceptance_ratio         = 0.2    # Target acceptance ration of new structures. Step_size will be updated to try and achieve this value;
const ar_buffer_zone           = 0.01   # Buffer zone limit around the target acceptance_ratio. Step_size will not be updated if it's currently inside the buffer zone;


# 8. Define the molecule secondary structure
const ss                       = "CEEEEEEECCCCHHHHHHHHHHHHCCCHHHHHHHHHCCCEEEEEEECHHHHHHHHHHHHHHCCEEEEC"


# 9. Define the contact/dihedral FBR settings
const contact_threshold        = 0.31
const contact_force_constant   = 10.0
const dihedral_force_constant  = 1e2
const dihedral_fb_width        = 10.0   # Width of the zone (in angles) where the restrictionm potential is zero (The dihedral angle can move half of this value each side)


# 10. Define the input files:
const input_pdb                = "../data/1ctf_no_sc.pdb"
const input_contact_map        = "../data/contact_map_raptorx.txt"
const input_amber_json         = "../data/1ctf_amber_top.json"


# 11. Define the output destination:
const log_destination          = stdout
const xyz_destination          = open("../out/trajectory.pdb", "w")
const out_destination          = open("../out/data.out", "w")
const outer_best_destination   = open("../out/outer_best_destination.pdb", "w")
const out_energy               = open("../out/energy.out", "w")
const status_destination       = [log_destination, out_destination]