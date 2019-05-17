# ----------------------------------------------------------------------------------------------------------
#                                                MACROS 
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
macro cbcall(cbs, state, driver_state)

    ex = quote
        for cb in $cbs
            if (getproperty(cb, :freq)>0) && (getproperty($driver_state, :step)%getproperty(cb, :freq)==0)
                getproperty(cb, :callback)($state, $driver_state)
            end
        end
    end
    return esc(ex)
end


@doc raw"""
    @Common.callback f::Function freq::Int64

(Macro) Create a [`CallbackObject`](@ref) with the given function `f` and output frequency `freq`. 

# Examples
```julia-repl
julia> @Common.callback 10 my_callback1
CallbackObject(freq=10, callback=my_callback1)

julia> @Common.callback my_callback1
CallbackObject(freq=1, callback=my_callback1)
```
"""
macro callback(args...)
    ex = :(Common.CallbackObject($(args...)))
    esc(ex)
end

macro copy(dst, src, components...)
    dst = esc(dst)
    src = esc(src)
    ex = quote; end
    for comp in components
        push!(ex.args, :(copy!($(dst).$(comp), $(src).$(comp))))
    end
    ex
end