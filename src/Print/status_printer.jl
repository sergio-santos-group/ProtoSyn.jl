@doc raw"""
    status(content::String, destination::Array{IO})

Print the `content` to all IO streams defined in `destination`.

# Examples
```julia-repl
julia> Print.status(@sprintf("%6d", step), log_file)

julia> Print.status(@sprintf("%6d", step), [log_file1, log_file2])
```
"""
function status(content::String, destination::Array{IO}, color::Union{Symbol,Int} = :normal)
    for io in destination
        printstyled(io, content, color=color)
        flush(io)
    end
end
function status(content::String, destination::IO = stdout)
    status(content, IO[destination])
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
    printstyled(destination, @sprintf("\n%20s: %10.3e\n", "⚡E_Total", energy.eTotal), color = :red)
    for (i, (comp, value)) in enumerate(energy.comp)
        write(destination, @sprintf "%2d: %15s: %10.3e\n" i comp value)
    end
    write(destination, "-"^31, "\n")
end