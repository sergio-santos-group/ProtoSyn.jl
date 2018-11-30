print_status_mc = @Common.callback print_sts_every_ic function cb_status(step::Int64, st::Common.State, dr::Drivers.MonteCarlo.Driver, ar::Float64, args...)
    Print.status(@sprintf("(%5s) %12d | ⚡E: %10.3e ▶️ ⚡Amber: %10.3e & ⚡Other: %10.3e | Dihedral: %8.2e | Crankshaft: %8.2e | AR: %.4f | T: %4.2f\n",
        "MC", step, st.energy.eTotal, st.energy.comp["amber"], st.energy.comp["other"], dihedral_mutator.step_size, crankshaft_mutator.step_size,
        ar, dr.temperature), status_destination)
end
print_status_sd = @Common.callback print_sts_every_min function cb_status(step::Int64, st::Common.State, dr::Drivers.SteepestDescent.Driver, max_force::Float64, gamma::Float64, args...)
    Print.status(@sprintf("(%5s) %12d | ⚡E: %10.3e ▶️ ⚡Amber: %10.3e & ⚡Other: %10.3e | Max Force: %10.3e | Gamma: %10.3e\n", "SD", step, st.energy.eTotal, st.energy.comp["amber"],
        st.energy.comp["other"], max_force, gamma), status_destination)
end

print_structure_ic = @Common.callback print_str_every_ic function cb_print_inner(step::Int64, st::Common.State, dr::Drivers.AbstractDriver, args...)
    Print.as_pdb(xyz_destination, st, metadata, step = step)
    flush(xyz_destination)
end
print_structure_min = @Common.callback print_str_every_min function cb_print_min(step::Int64, st::Common.State, dr::Drivers.AbstractDriver, args...)
    Print.as_pdb(xyz_destination, st, metadata, step = step)
    flush(xyz_destination)
end

outer_best_energy = Inf
print_outer_best = @Common.callback 1 function cb_print_outer(step::Int64, st::Common.State, dr::Drivers.AbstractDriver, args...)
    global outer_best_energy
    if st.energy.eTotal < outer_best_energy
        Print.as_pdb(outer_best_destination, st, metadata, step = step)
        flush(outer_best_destination)
        println(@sprintf "(ILSRR) New outer best defined: ⚡E: %10.3e (old) ▶️ %10.3e (new)" outer_best_energy st.energy.eTotal)
        outer_best_energy = st.energy.eTotal
    end
end

adjust_step_size = @Common.callback 1 function cb_adjust_step_size(step::Int64, st::Common.State, dr::Drivers.MonteCarlo.Driver, ac, args...)
    d::Float64 = 1.0
    if ac > acceptance_ratio + ar_buffer_zone
        d = 1.0005
    elseif ac < acceptance_ratio - ar_buffer_zone
        d = 0.9995
    end
    dihedral_mutator.step_size = max(min_step_size, min(dihedral_mutator.step_size * d, π))
    crankshaft_mutator.step_size = max(min_step_size, min(crankshaft_mutator.step_size * d, π))
end

quench_temperature = @Common.callback 1 function cb_adjust_temperature(step::Int64, st::Common.State, dr::Drivers.MonteCarlo.Driver, args...)
    dr.temperature = -(temperature/(dr.n_steps - 1))*step + temperature
end

reset_temperature = @Common.callback 1 function cb_reset_temperature(step::Int64, st::Common.State, dr::Drivers.ILSRR.Driver, args...)
    dr.inner_cycle_driver.temperature = temperature
end

reset_step_size = @Common.callback 1 function cb_reset_step_size(step::Int64, st::Common.State, dr::Drivers.ILSRR.Driver, args...)
    dihedral_mutator.step_size = dihedral_step_size
    crankshaft_mutator.step_size = crankshaft_step_size
end