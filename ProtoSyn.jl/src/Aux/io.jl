@doc raw"""
    read_JSON(i_file::String)

Read a JSON file and return a dictionary with the content.

# Examples
```julia-repl
julia> Aux.read_JSON("i_file.json")
Dict{Any, Any}()
```
"""
function read_JSON(i_file::String)::Dict

    open(i_file, "r") do f
        json_txt = read(f, String)
        return JSON.parse(json_txt)
    end
end