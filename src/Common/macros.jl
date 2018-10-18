# ----------------------------------------------------------------------------------------------------------
#                                                CALLBACK 
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

    ex = quote
        for cb in $cbs
            if (getproperty(cb, :freq)>0) && ($step%getproperty(cb, :freq)==0)
                getproperty(cb, :callback)($step, $(args...))
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
julia> @Common.callback my_callback1 10
CallbackObject(callback=my_callback1, freq=10)
```
"""
macro callback(f, freq)
    ex = :(CallbackObject($f, $freq))
    esc(ex)
end

# ----------------------------------------------------------------------------------------------------------
#                                                Energy
@doc raw"""
    @Common.faggregator name::String f::function Vararg::Any

(Macro) Aggregate multiple functions `f` in a single variable `name`. `Vararg` contains the arguments used by function `f`

# Examples
```julia-repl
julia> @faggregator myeval f top1
@faggregator myeval g top2
@faggregator myeval h top3

energy = myevalf(state, false)
```
"""
macro faggregator(name, f, args...)
    cachename = Symbol("##", "_faggregator_cache")
    mod = @__MODULE__
    
    cache = isdefined(mod, cachename) ?
             getfield(mod, cachename) :
             Core.eval(mod, :(const $cachename = Dict{Symbol, Expr}()))

    if haskey(cache, name, )
        ex = :(energy += $f($(args...), s, b))
        expr = cache[name]
        pos = length(expr.args[2].args[2].args) - 2
        insert!(expr.args[2].args[2].args, pos, ex)
    else
        cache[name] = quote
            function $name(s::State, b::Bool)
                energy = $f($(args...), s, b)
                s.eTotal = energy
                return energy
            end
        end
    end
    esc(cache[name])
end