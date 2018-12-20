# PRINT STATUS ----
print_status_mc = @Common.callback print_sts_every_ic function _print_status_mc(step::Int64, st::Common.State, dr::Drivers.MonteCarlo.Driver, ar::Float64, args...)
    Print.status(@sprintf("(%5s) %12d | ⚡E: %10.3e ▶️ ⚡Amber: %10.3e & ⚡Other: %10.3e | Dihedral: %8.2e | Crankshaft: %8.2e | AR: %.4f | T: %4.2f\n",
        "MC", step, st.energy.eTotal, st.energy.comp["amber"], st.energy.comp["other"], reg_dh_mutator.step_size, reg_cs_mutator.step_size,
        ar, dr.temperature), status_destination)
end

print_status_refnm = @Common.callback print_sts_every_rfmm function _print_status_refnm(step::Int64, st::Common.State, dr::Drivers.MonteCarlo.Driver, ar::Float64, args...)
    Print.status(@sprintf("(%5s) %12d | ⚡E: %10.3e ▶️ ⚡Amber: %10.3e & ⚡Other: %10.3e | DH: %8.2e | CS: %8.2e | BR: %8.2f | AR: %.4f | T: %4.2f\n",
        "REFNM", step, st.energy.eTotal, st.energy.comp["amber"], st.energy.comp["other"], soft_dh_mutator.step_size, soft_cs_mutator.step_size, soft_br_mutator.step_size,
        ar, dr.temperature), status_destination, :green)
end

print_status_init = @Common.callback print_sts_every_min function _print_status_init(step::Int64, st::Common.State, dr::Drivers.SteepestDescent.Driver, max_force::Float64, gamma::Float64, args...)
    Print.status(@sprintf("(%5s) %12d | ⚡E: %10.3e ▶️ ⚡Amber: %10.3e & ⚡Other: %10.3e | Max Force: %10.3e | Gamma: %10.3e\n", "I_SD", step, st.energy.eTotal,
        st.energy.comp["amber"], st.energy.comp["other"], max_force, gamma), status_destination, 11)
end

print_status_loop = @Common.callback print_sts_every_min function _print_status_loop(step::Int64, st::Common.State, dr::Drivers.SteepestDescent.Driver, max_force::Float64, gamma::Float64, args...)
    Print.status(@sprintf("⤷(%4s) %12d | ⚡E: %10.3e ▶️ ⚡Amber: %10.3e & ⚡Other: %10.3e | Max Force: %10.3e | Gamma: %10.3e\n", "LOOP", step, st.energy.eTotal,
        st.energy.comp["amber"], st.energy.comp["other"], max_force, gamma), status_destination)
end

print_status_smpl = @Common.callback print_sts_every_min function _print_status_smpl(step::Int64, st::Common.State, dr::Drivers.SteepestDescent.Driver, max_force::Float64, gamma::Float64, args...)
    Print.status(@sprintf("⤷(%4s) %12d | ⚡E: %10.3e ▶️ ⚡Amber: %10.3e & ⚡Other: %10.3e | Max Force: %10.3e | Gamma: %10.3e\n", "SMPL", step, st.energy.eTotal,
        st.energy.comp["amber"], st.energy.comp["other"], max_force, gamma), status_destination)
end

# PRINT STRUCTURE ----
print_structure = @Common.callback print_str_every_oc function _print_structure(step::Int64, st::Common.State, dr::Drivers.AbstractDriver, args...)
    Print.as_pdb(xyz_destination, st, metadata, step = step)
    flush(xyz_destination)
end

print_structure_rfnm = @Common.callback print_str_every_rfnm function _print_structure_rfnm(step::Int64, st::Common.State, dr::Drivers.AbstractDriver, args...)
    Print.as_pdb(xyz_destination, st, metadata, step = step)
    flush(xyz_destination)
end


#PRINT ENERGY
print_energy = @Common.callback print_ene_every_ic function _print_energy(step::Int64, st::Common.State, dr::Drivers.AbstractDriver, args...)
    println(out_ene, st.energy.eTotal)
    flush(out_ene)
end

print_energy_components = @Common.callback print_ene_every_ic function _print_energy_components(step::Int64, st::Common.State, dr::Drivers.AbstractDriver, args...)
    
    title = ""
    try
        title = args[2]
    catch
        nothing
    end
    write(out_ene, @sprintf("\n%6s %6s ", title, step))
    for component in energy_components_list
        if component == "eTotal"
            write(out_ene, @sprintf("%11.4e ", st.energy.eTotal))
        elseif component in keys(st.energy.comp)
            write(out_ene, @sprintf("%11.4e ", st.energy.comp[component]))
        else
            write(out_ene, @sprintf("%11s ", "NaN"))
        end
    end
    flush(out_ene)
end


outer_best_energy = Inf
print_structure_best = @Common.callback 1 function _print_structure_best(step::Int64, st::Common.State, dr::Drivers.AbstractDriver, args...)
    global outer_best_energy
    if st.energy.eTotal < outer_best_energy
        Print.as_pdb(outer_best_destination, st, metadata, step = step)
        flush(outer_best_destination)
        printstyled(@sprintf("(ILSRR) New outer best defined: ⚡E: %10.3e (old) ▶️ %10.3e (new)\n", outer_best_energy, st.energy.eTotal), color = :red)
        outer_best_energy = st.energy.eTotal
        
        #Print energy
        _print_energy_components(step, st, dr, 0.0, "(BEST)")
    end
end

# SIMULATION ADJUSTMENTS ---
adjust_step_size = @Common.callback 1 function _adjust_step_size(step::Int64, st::Common.State, dr::Drivers.MonteCarlo.Driver, ac::Float64, args...)
    d::Float64 = 1.0
    if ac > acceptance_ratio + ar_buffer_zone
        d = 1.0 + step_size_adjust_scale
    elseif ac < acceptance_ratio - ar_buffer_zone
        d = 1.0 - step_size_adjust_scale
    end
    reg_dh_mutator.step_size = max(min_step_s, min(reg_dh_mutator.step_size * d, π))
    reg_cs_mutator.step_size = max(min_step_s, min(reg_cs_mutator.step_size * d, π))
end

adjust_step_size_rfnm = @Common.callback 1 function _adjust_step_size_refnm(step::Int64, st::Common.State, dr::Drivers.MonteCarlo.Driver, ac::Float64, args...)
    d::Float64 = 1.0
    if ac > acceptance_ratio + ar_buffer_zone
        d = 1.05
    elseif ac < acceptance_ratio - ar_buffer_zone
        d = 0.95
    end
    # soft_dh_mutator.step_size = max(min_step_s, min(soft_dh_mutator.step_size * d, π))
    # soft_cs_mutator.step_size = max(min_step_s, min(soft_cs_mutator.step_size * d, π))
    soft_br_mutator.step_size = max(min_step_s, min(soft_br_mutator.step_size * d, π))
end

quench_temperature = @Common.callback 1 function cb_adjust_temperature(step::Int64, st::Common.State, dr::Drivers.MonteCarlo.Driver, args...)
    dr.temperature = - (i_search_temp_init / (dr.n_steps - 1)) * step + i_search_temp_init
end

reset_temperature = @Common.callback 1 function cb_reset_temperature(step::Int64, st::Common.State, dr::Drivers.ILSRR.Driver, args...)
    dr.inner_cycle_driver.temperature = i_search_temp_init
end

reset_step_size = @Common.callback 1 function cb_reset_step_size(step::Int64, st::Common.State, dr::Drivers.ILSRR.Driver, args...)
    reg_dh_mutator.step_size = reg_dh_step_s
    reg_cs_mutator.step_size = reg_cs_step_s
end