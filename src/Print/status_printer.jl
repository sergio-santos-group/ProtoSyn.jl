@doc raw"""
    status(content::String, destination::Array{IO})

Print the `content` to all IO streams defined in `destination`.

# Examples
```julia-repl
julia> Print.status(@sprintf("%6d", step), log_file)

julia> Print.status(@sprintf("%6d", step), [log_file1, log_file2])
```
"""
function status(content::String, destination::Array{IO})
    for io in destination
        write(io, content)
        flush(io)
    end
end
function status(content::String, destination::IO)
    status(content, [destination])
end