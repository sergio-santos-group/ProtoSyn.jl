using ProtoSyn
using Printf
using LinearAlgebra

include("ilsrr_config.jl")

state, metadata        = Common.load_from_pdb(input_pdb)
amber_topology         = Forcefield.Amber.load_from_json(input_amber_json)
Common.fix_proline!(state, metadata.dihedrals)
Common.apply_ss!(state, metadata, ss)
contact_restraints     = Forcefield.Restraints.load_distance_restraints_from_file(input_contact_map, metadata, k = contact_force_constant, threshold = contact_threshold)
dihedral_restraints    = Forcefield.Restraints.lock_block_bb(metadata, k = dihedral_force_constant)

nb_dihedrals           = filter(x -> x.residue.ss == Common.SS.COIL, metadata.dihedrals)
nb_phi_dihedrals       = filter(x -> x.dtype == Common.DIHEDRAL.phi, nb_dihedrals) # For crankshaft movements
dihedral_mutator       = Mutators.Dihedral.DihedralMutator(nb_dihedrals, () -> (randn() * dihedral_mutator.step_size), dihedral_p_mut, dihedral_step_size)
crankshaft_mutator     = Mutators.Crankshaft.CrankshaftMutator(nb_phi_dihedrals, () -> (randn() * crankshaft_mutator.step_size), crankshaft_p_mut, crankshaft_step_size)
dihedral_perturbator   = Mutators.Dihedral.DihedralMutator(nb_dihedrals, () -> (randn() * dihedral_mutator.step_size), dihedral_p_mut_prt, dihedral_step_size_prt)
crankshaft_perturbator = Mutators.Crankshaft.CrankshaftMutator(nb_phi_dihedrals, () -> (randn() * crankshaft_mutator.step_size), crankshaft_p_mut_prt, crankshaft_step_size_prt)

include("ilsrr_callbacks.jl")

function my_evaluator!(st::Common.State, do_forces::Bool)
    fill!(st.forces, zero(Float64))
    # other_energy  = Forcefield.Restraints.evaluate!(contact_restraints, st, do_forces=do_forces)
    other_energy = Forcefield.Restraints.evaluate!(dihedral_restraints, st, do_forces=do_forces)
    amber_energy = Forcefield.Amber.evaluate!(amber_topology, st, cut_off=1.2, do_forces=do_forces)
    # amber_energy = 0
    # other_energy = 0
    st.energy.comp["other"] = other_energy
    st.energy.eTotal = amber_energy + other_energy
    println(@sprintf "Trying: %10.3e" st.energy.eTotal)
    return amber_energy + other_energy
end

sampler_sd_driver = Drivers.SteepestDescent.Driver(my_evaluator!, min_n_steps)
initial_sd_driver = Drivers.SteepestDescent.Driver(my_evaluator!, init_min_n_steps, print_status_sd)
function my_sampler!(st::Common.State)
    Mutators.Dihedral.run!(st, dihedral_mutator)
    Mutators.Crankshaft.run!(st, crankshaft_mutator)
    # sampler_sd_driver.run!(st, sampler_sd_driver)
    Print.as_pdb(xyz_destination, state, metadata)
end

function my_pertubator!(st::Common.State)
    dm = Mutators.Dihedral.run!(st, dihedral_perturbator)
    cm = Mutators.Crankshaft.run!(st, crankshaft_perturbator)
    println("(ILSRR) Performing structure perturbation ▶️ Dihedral: $dm | Crankshaft: $cm")
    sampler_sd_driver.run!(st, initial_sd_driver)
end

inner_cycle_driver = Drivers.MonteCarlo.Driver(my_sampler!, my_evaluator!, temperature, n_inner_steps, print_status_mc, quench_temperature, print_structure_ic, adjust_step_size)
ilsrr_driver = Drivers.ILSRR.Driver(inner_cycle_driver, my_evaluator!, my_pertubator!, temperature, n_outer_steps, reset_temperature, reset_step_size, print_outer_best)
# initial_sd_driver.run!(state, initial_sd_driver, print_structure_min)
# ilsrr_driver.run!(state, ilsrr_driver)



function calculate_dihedral(a1::Vector{Float64}, a2::Vector{Float64}, a3::Vector{Float64}, a4::Vector{Float64})::Float64

    v12 = a2 - a1
    v23 = a3 - a2
    v34 = a4 - a3
    v123 = cross(v12, v23)
    v234 = cross(v23, v34)
    return atan(dot(cross(v123, v234), v23)/sqrt(dot(v23, v23)), dot(v123, v234))
end

function get_FBR_from_dihedral(dihedral::Common.Dihedral, restraint_list::Vector{Forcefield.Restraints.DihedralFBR})::Union{Nothing, Forcefield.Restraints.DihedralFBR}

    for restraint in restraint_list
        if (restraint.a1 == dihedral.a1) && (restraint.a2 == dihedral.a2) && (restraint.a3 == dihedral.a3) && (restraint.a4 == dihedral.a4)
            return restraint
        end
    end
    return nothing
end

for dihedral in metadata.dihedrals
    fbr = get_FBR_from_dihedral(dihedral, dihedral_restraints)
    if dihedral.dtype < Common.DIHEDRAL.omega && fbr != nothing
        angle = calculate_dihedral(state.xyz[dihedral.a1, :], state.xyz[dihedral.a2, :], state.xyz[dihedral.a3, :], state.xyz[dihedral.a4, :])
        println(@sprintf "%s(%3d-%3d-%3d-%3d)%s -> %7.2f (%7.2f|%7.2f|%7.2f|%7.2f)" dihedral.dtype dihedral.a1 dihedral.a2 dihedral.a3 dihedral.a4 dihedral.residue.name rad2deg(angle) rad2deg(fbr.r1) rad2deg(fbr.r2) rad2deg(fbr.r3) rad2deg(fbr.r4))
        # println(" $(rag2deg(fbr.r1))-$(rag2deg(fbr.r2))-$(rag2deg(fbr.r3))-$(rag2deg(fbr.r4))")
    end
end

# for dihedral in dihedral_restraints
#     println(dihedral)
# end