# ----------------------------------------------------------------------------------------------------------
#                                                CALLBACK 

@doc raw"""
    CallbackObject(freq::Int64, callback::Function)

Define the callback function parameters.

# Arguments
- `freq`: Frequency (in steps) that the callback function is called.
- `callback`: Actual callback function.

# Examples
```julia-repl
julia> Common.CallbackObject(1, Print.as_xyz)
CallbackObject(freq=1, callback=Print.as_xyz)

See also: [`Print.as_xyz`](@ref Print) [`cbcall`](@ref)
```
"""
mutable struct CallbackObject
    freq::Int64
    callback::Function
end
Base.show(io::IO, b::CallbackObject) = print(io, "CallbackObject(freq=$(b.freq), callback=$(string(b.callback)))")