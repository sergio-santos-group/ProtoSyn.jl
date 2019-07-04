# ------------------ Aux functions -------------------------

function print_energy_components(driver::String, step::Int64, energy::Common.Energy, temp::Float64)::String
    components = [:amber, :contacts, :sol, :hb]
    s = @sprintf("%10s %6d %11.4e", driver, step, energy.total)
    for component in components
        if component in keys(energy.components)
            s = join([s, @sprintf("%11.4e", energy.components[component])], "  ")
        else
            s = join([s, @sprintf("%11s", "----------")], "  ")
        end
    end
    s = join([s, @sprintf("| %6.3f", temp)], "  ")
    return s
end

function print_energy_keys(driver::String, step::Int64)::String
    components = [:amber, :contacts, :sol, :hb]
    s = @sprintf("%10s %6d %11s", driver, step, "TOTAL")
    for key in components
        s = join([s, @sprintf("%11s", uppercase(string(key)))], "  ")
    end
    s = join([s, "  | TEMPERATURE"])
    return s
end