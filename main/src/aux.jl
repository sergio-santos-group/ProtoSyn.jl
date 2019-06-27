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

function printcredits()
    printstyled("\n    P R O T O S Y N   M A I N    A L G O R I T H M\n", color = :bold)
    printstyled("  > Created by Sérgio M. Santos & José M. Pereira\n", color = :blue)
    printstyled("  > Universidade de Aveiro, Portugal, 2019\n", color = :blue)
    printstyled("  > The ProtoSyn.jl package is licensed under the MIT 'Expat'\n", color = :blue)
    printstyled("  > Copyright (c) 2018: sergio-santos-group\n", color = :blue)
    printstyled("  > Requires ProtoSyn library v0.2 or higher\n", color = :blue)
    # printstyled("  > If useful, please cite [ARTICLE]\n", color = :blue)
    # printstyled("  > For more information contact [CONTACT]\n", color = :blue)
    println("")
end