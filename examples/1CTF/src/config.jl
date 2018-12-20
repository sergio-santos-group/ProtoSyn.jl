# Algorithm length (number of steps):
const n_search_refine_cycles   = 10000  # Number of cycles between Random Search / Refinement stages;
const n_i_search_steps         = 6000   # Number of steps on each Monte Carlo Search inner cycle;
const n_o_search_steps         = 100    # Number of ILSRR outer cycles to perform;
const n_refine_steps           = 10000  # Number of steps for Monte Carlo Refinement;

const n_init_min_steps         = 10000  # Number of steps for Steepest Descent minimization (initial long minimization);
const n_refine_min_steps       = 500    # Number of steps for Steepest Descent minimization (every step of Monte Carlo Refinement);
const n_loop_closer_steps      = 500    # Number of steps for Steepest Descent in loop closure (After Blockrot movements in Monte Carlo Refinement);


# Output frequency:
const print_sts_every_ic       = 1000   # Print system status every n steps (Inner cycle);
const print_sts_every_min      = 10     # Print status every n steps in Steepest Descent minimization;
const print_sts_every_rfmm     = 1      # Print status every n steps in Monte Carlo Refinement;

const print_str_every_oc       = 1      # Print structure to files every n steps;
const print_str_every_rfnm     = 0      # Print structure to files every n steps in Monte Carlo Refinement;

const print_ene_every_ic       = 1000   # Print energy information every n steps;

# Slope anaylis parameters
const evaluate_slope_every     = 500    # Calculate the energy decrease slope of the last n steps;
const evaluate_slope_threshold = -1.0   # Convergence criteria: if the energy decrease slope is smaller than this threshold, the system moves on;
const continue_after_n_attemps = 15     # Exit the ISLRR algorithm after n jumps failed to produce a new homebase; 

# Mutators step size:
const min_step_s               = 0.01   # Minimum step_size allowed;
const soft_dh_step_s           = π/32   # Step size for dihedral movements in the Monte Carlo Refinement;
reg_dh_step_s                  = π/8    # (Initial) Step size for dihedral movements in the Monte Carlo Search (ILSRR inner cycle);
const hard_dh_step_s           = π/1    # Step size for dihedral movements in the ILSRR perturbator;

const soft_cs_step_s           = π/32   # Step size for crankshaft movements in the Monte Carlo Refinement;
reg_cs_step_s                  = π/8    # (Initial) Step size for crankshaft movements in the Monte Carlo Search (ILSRR inner cycle);
const hard_cs_step_s           = π/1    # Step size for crankshaft movements in the ILSRR perturbator;

const soft_br_step_s           = π/2    # Step size for blockrot movements in the Monte Carlo Refinement;
const soft_br_tr_step_s        = 0.2    # Step size for translation in blockrot movements;

# Mutators probabilities:
const soft_dh_p_mut            = 0.025  # Mutation probability for dihedral movements in the Monte Carlo Refinement;
const reg_dh_p_mut             = 0.025  # Mutation probability for dihedral movements in the Monte Carlo Search (ILSRR inner cycle);
const hard_dh_p_mut            = 0.1    # Mutation probability for dihedral movements in the ILSRR perturbator;

const soft_cs_p_mut            = 0.007  # Mutation probability for crankshaft movements in the Monte Carlo Refinement;
const reg_cs_p_mut             = 0.007  # Mutation probability for crankshaft movements in the Monte Carlo Search (ILSRR inner cycle);
const hard_cs_p_mut            = 0.1    # Mutation probability for crankshaft movements in the ILSRR perturbator;

const soft_br_p_mut            = 0.15   # Mutation probability for blockrot movements in the Monte Carlo Refinement;
const soft_br_tries            = 10     # Number of attempts to rotate a block after it was accepted by p_mut;

# Steepest Descent parameters (convergence criteria and max step size):
const max_step                 = 1e-2   # Maximum step allowed in Steepest Descent;
const f_tol                    = 0.1    # Minimum force tolerated. The minimization is considered converged when forces are below the defined value;

# Forcefield and energy function parameters:
const nonbonded_cut_off        = 1.2    # Cutoff (in nm) for nonbonded interactions calculation;
const λ_eContact               = 1.0   # Scalar constant for contact energy;
const λ_eContact_rfnm          = 2000.0 # Scalar constant for contact energy in the Refinement stage;
const λ_eSol                   = 1.0   # Scalar constant for solvation energy;
const λ_eSol_rfnm              = 10.0   # Scalar constant for solvation energy in the Refinement stage;

# System temperature(s):
i_search_temp_init             = 7.0    # Temperature for Metropolis acceptance in inner cycles (0.0 performs a Random Search where only better structures are accepted);
const o_search_temp_static     = 10.0   # Temperature for ILSRR acceptance in outer cycles (0.0 performs a Random Search where only better structures are accepted);
const refine_temp_static       = 1.0    # Temperature for Metropolois acceptance in refinement (0.0 performs a Random Search where only better structures are accepted);

# Target acceptance ratio for Monte Carlo Drivers and respective buffer zone:
const acceptance_ratio         = 0.24   # Target acceptance ration of new structures. Step_size will be updated to try and achieve this value;
const ar_buffer_zone           = 0.002  # Buffer zone limit around the target acceptance_ratio. Step_size will not be updated if it's currently inside the buffer zone;
const step_size_adjust_scale   = 5e-5   # Amount of scaling when adjust step size (in both directions)

# Contact FBR settings
const contact_threshold        = 0.00   # Consider only contact with certainty percentages above the defined threshold;
const contact_force_constant   = 2000.0 # Force constant for contact distance FBRs;
const contact_min_distance     = 0.8    # Minimum distance (in nm) between the atoms involved in the contact;

# Dihedral FBR settings
const dihedral_force_constant  = 10.0  # Force constant for dihedral FBRs; 
const dihedral_fb_width        = 10.0   # Width of the zone (in angles) where the restrictionm potential is zero (The dihedral angle can move half of this value each side)


# Secondary structure (C - Coil | H - α-Helix | E - β-sheet)
const ss                       = "CEEEEEEECCCCHHHHHHHHHHHHCCCHHHHHHHHHCCCEEEEEEECHHHHHHHHHHHHHHCCEEEEC"

# Input files:
const input_pdb                = "data/1ctf_no_sc.pdb"
const input_contact_map        = "data/contact_map_raptorx.txt"
const input_amber_json         = "data/1ctf_amber_top_no_sc.json"
const ref_native               = "data/1ctf_native.pdb"
const conf_pdb                 = "sample.pdb"


# Output destination:
const log_destination          = stdout
const xyz_destination          = open("out/trajectory.pdb", "w")
const out_destination          = open("out/data.out", "w")
const outer_best_destination   = open("out/outer_best_destination.pdb", "w")
const status_destination       = [log_destination, out_destination]
const out_ene                  = open("out/energy.out", "w")
const energy_components_list   = ["eTotal", "eBond", "eAngle", "eDihedral", "eCoulomb", "eCoulomb14", "eLJ", "eLJ14", "eDistanceFBR", "eSol", "eDihedralFBR"]
println(out_ene, @sprintf("%6s %6s %11s %11s %11s %11s %11s %11s %11s %11s %11s %11s %11s", "", "Step", "eTotal", "eBond", "eAngle", "eDihedral", "eCoulomb", "eCoulomb14", "eLJ", "eLJ14", "eContacts", "eSol", "eDihedralFBR"))