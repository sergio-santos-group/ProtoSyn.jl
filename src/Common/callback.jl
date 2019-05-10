# ----------------------------------------------------------------------------------------------------------
#                                                CALLBACK 

@doc raw"""
    CallbackObject(freq::Int64, callback::Function)

Define the callback function parameters.

# Arguments
- `freq`: Frequency (in steps) that the callback function is called.
- `callback`: Actual callback function. This function should have the following signature:
```
callback(step::Int64, st::Common.State, dr::Drivers.MonteCarlo.Driver, args...)
```

# Examples
```julia-repl
julia> Common.CallbackObject(100, Print.as_xyz)
CallbackObject(freq=100, callback=Print.as_xyz)

julia> Common.CallbackObject(Print.as_xyz)
CallbackObject(freq=1, callback=Print.as_xyz)
```
See also: [`Print.as_xyz`](@ref Print) [`@cbcall`](@ref)
"""
Base.@kwdef mutable struct CallbackObject{F<:Function} <: Abstract.CallbackObject
    freq::Int = 1
    callback::F
end
Base.show(io::IO, b::CallbackObject) = begin
    print(io, "CallbackObject(freq=$(b.freq), callback=$(string(b.callback)))")
end