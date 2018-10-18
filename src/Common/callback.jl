# ----------------------------------------------------------------------------------------------------------
#                                                CALLBACK 

@doc raw"""
    CallbackObject(freq::Int64, callback::Function)

Define the callback function parameters.

# Arguments
- `freq`: Frequency (in steps) that the callback function is called.
- `callback`: Actual callback function. This function should have the following signature:
```
callback(step::Int64, st::Common.State, dr::Drivers.MonteCarlo.MonteCarloDriver, args...)
```

# Examples
```julia-repl
julia> Common.CallbackObject(100, Print.as_xyz)
CallbackObject(callback=Print.as_xyz, freq=100)

julia> Common.CallbackObject(Print.as_xyz)
CallbackObject(callback=Print.as_xyz, freq=1)
```
See also: [`Print.as_xyz`](@ref Print) [`@cbcall`](@ref)
"""
mutable struct CallbackObject
    callback::Function
    freq::Int64
end
CallbackObject(callback::Function; freq::Int64 = 1) = CallbackObject(callback, freq)
Base.show(io::IO, b::CallbackObject) = print(io, "CallbackObject(callback=$(string(b.callback)), freq=$(b.freq))")