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

See also: [`Print.as_xyz`](@ref Print)
```
"""
mutable struct CallbackObject
    freq::Int64
    callback::Function
end
Base.show(io::IO, b::CallbackObject) = print(io, "CallbackObject(freq=$(b.freq), callback=$(string(b.callback)))")

@doc raw"""
    @Common.cbcall callbacks::Tuple{CallbackObject, N} step::Int64 Vararg::Any

(Macro) Call the CallbackObject.function of each [`CallbackObject`](@ref) in the `callbacks` Tuple{CallbackObject, N} independently.
Each CallbackObject.function is ran depending on the defined CallbackObject.freq and the given `step`.
Vararg holds all the arguments necessary to run the callback function itself.

# Examples
```julia-repl
julia> @Common.cbcall (callback_object1, callback_object2) 1
```
"""
macro cbcall(cbs, step, args...)
    tmp = eval(cbs)

    if length(tmp) == 0
        ex = :nothing
    else
        ex = quote
            for cb in $cbs
                if (getproperty(cb, :freq)>0) && ($step%getproperty(cb, :freq)==0)
                    getproperty(cb, :callback)($step, $(args...))
                end
            end
        end
    end
    return esc(ex)
end