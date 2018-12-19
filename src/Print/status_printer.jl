@doc raw"""
    status(content::String[, destination::Array{IO} = stdout, color::Union{Symbol, Int} = :normal])

Print the `content` to all IO streams defined in `destination`, in the defined `color`.

# Examples
```julia-repl
julia> Print.status(@sprintf("%6d", step), log_file, :red)

julia> Print.status(@sprintf("%6d", step), [log_file1, log_file2])
```
"""
function status(content::String, destination::Array{IO}, color::Union{Symbol, Int} = :normal)
    for io in destination
        printstyled(io, content, color=color)
        flush(io)
    end
end
function status(content::String; destination::IO = stdout, color::Union{Symbol, Int} = :normal)
    status(content, IO[destination], color)
end


@doc raw"""
    energy_by_component(energy::Common.Energy,[ destination::IO = stdout])

Print the `energy` content to the IO stream defined in `destination`, by component.

# Examples
```julia-repl
julia> Print.energy_by_component(state.energy, file_name)
```
"""
function energy_by_component(energy::Common.Energy, destination::IO = stdout)

    write(destination, "-"^31)
    printstyled(destination, @sprintf("\n%20s: %10.3e\n", "⚡E_Total", energy.eTotal), color = :blue)
    for (i, (comp, value)) in enumerate(energy.comp)
        write(destination, @sprintf "%2d: %15s: %10.3e\n" i comp value)
    end
    write(destination, "-"^31, "\n")
end

# TODO: Add documentation
function energy_by_component(ref_energy::Common.Energy, sample_energy::Common.Energy; ref_title::String = "", sample_title = "", destination::IO = stdout)

    write(destination, "-"^55)
    write(destination, @sprintf("\n%31s %10s", ref_title, sample_title))
    printstyled(destination, @sprintf("\n%20s: %10.3e %10.3e ", "⚡E_Total", ref_energy.eTotal, sample_energy.eTotal), color = :blue)
    diff = ref_energy.eTotal - sample_energy.eTotal
    c = diff >= 0 ? :red : :green
    printstyled(destination, @sprintf("(%+10.3e)\n", diff), color = c)
    for (i, (comp, value)) in enumerate(ref_energy.comp)
        write(destination, @sprintf "%2d: %15s: %10.3e " i comp value)
        write(destination, @sprintf "%10.3e " sample_energy.comp[comp])
        diff = value - sample_energy.comp[comp]
        c = diff >= 0 ? :red : :green
        printstyled(destination, @sprintf("(%+10.3e)\n", diff), color = c)
    end
    write(destination, "-"^55, "\n")
end