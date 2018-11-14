# ----------------------
# DEFINE THE CALLBACKS:
# ----------------------

# 1. Print status ------------------------------------------------------------------------------------ 
dr_type = Union{Drivers.MonteCarlo.MonteCarloDriver, Drivers.SteepestDescent.SteepestDescentDriver}
function cb_status(step::Int64, st::Common.State, dr::Drivers.MonteCarlo.MonteCarloDriver, ar::Float64, mov_count::Dict{String, Int64}, args...)
    Print.status(@sprintf("(%5s) %12d | ⚡E: %10.3e ▶️ ⚡Amber: %10.3e & ⚡Other: %10.3e | Dihedral: %2d ▶️ %8.2e | Crankshaft: %2d ▶️ %8.2e | AR: %.4f | T: %4.2f\n",
        "MC",
        step, st.energy.eTotal,
        st.energy.comp["amber"],
        st.energy.comp["eContact"],
        mov_count["d"],
        dihedral_mutator.step_size,
        mov_count["c"],
        crankshaft_mutator.step_size,
        ar,
        dr.temperature),
        status_destination)
end
print_status_mc = Common.CallbackObject(print_sts_every_ic, cb_status) 
print_status_sd = @Common.callback print_sts_every_min function cb_status(step::Int64, st::Common.State, dr::Drivers.SteepestDescent.SteepestDescentDriver, args...)
    Print.status(@sprintf("(%5s) %12d | ⚡E: %10.3e | Max Force: %10.3e | Gamma: %10.3e\n", "SD", step, st.energy.eTotal, args[1], args[2]),
        status_destination)
end


# 2. Check if the structure is better than the current inner_best and save it -------------------------
save_inner_best = @Common.callback 1 function cb_save(step::Int64, st::Common.State, dr::dr_type, args...)
    global inner_best
    if st.energy.eTotal < inner_best.energy.eTotal
        inner_best = deepcopy(st)
        Print.as_pdb(best_destination, st, metadata, title = "$step Inner best", step = step)
    end
end


# 3. Adjust the step_size so that, on average, the acceptance_ration is as defined initially -----------
adjust_step_size = @Common.callback updt_step_size_every function cb_adjust_step_size(step::Int64, st::Common.State, dr::Drivers.MonteCarlo.MonteCarloDriver, ac, args...)
    d::Float64 = 1.0
    if ac > acceptance_ratio + ar_buffer_zone
        d = 1.0005
    elseif ac < acceptance_ratio - ar_buffer_zone
        d = 0.9995
    end
    dihedral_mutator.step_size = max(min_step_size, min(dihedral_mutator.step_size * d, π))
    crankshaft_mutator.step_size = max(min_step_size, min(crankshaft_mutator.step_size * d, π))
end


# 4. Print current structure to a PDB file -------------------------------------------------------------
print_structure = @Common.callback print_str_every_ic function cb_print(step::Int64, st::Common.State, dr::dr_type, args...)
    Print.as_pdb(xyz_destination, st, metadata, step = step)
end


# 5. Adjust temperature according to a given annealing function
adjust_temperature = @Common.callback 1 function cb_adjust_temperature(step::Int64, st::Common.State, dr::Drivers.MonteCarlo.MonteCarloDriver, args...)
    dr.temperature = -(temperature/n_inner_steps)*step + temperature
end