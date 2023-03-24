using Base.Cartesian

export @reduce

"""
    @reduce(n::Int, op::Symbol, ex::Expr)

        Wrapping macro for @ncall. Is equivalent and generates
        op(ex_1, ..., ex_n)

    Example:

        @reduce 2 (+) u -> v_u^2
        > Generates +(v_1^2, v_2^2)
"""
macro reduce(n::Int, op::Symbol, ex::Expr)
    esc(:(@ncall($n, $op, $ex)))
end

@doc """
> @macroexpand @dot u a_u b_u
:(a_1 * b_1 + a_2 * b_2 + a_3 * b_3)
"""
macro dot(s::Symbol, v1::Symbol, v2::Symbol)
    esc(:(@ncall(3,+,$s -> $v1*$v2)))
end



#macro d(p::Symbol, n::Int)
#    s = Symbol(p, "_", (n-1)%3+1)
#    esc(:($s))
#end

cycle(n::Int) = (n-1)%3+1


@doc """
> @macroexpand @cross u out_u a_u b_u
quote
    out_1 = a_2 * b_3 - a_3 * b_2
    out_2 = a_3 * b_1 - a_1 * b_3
    out_3 = a_1 * b_2 - a_2 * b_1
end
"""
macro cross(s::Symbol, v::Union{Symbol, Expr}, v1::Union{Symbol, Expr}, v2::Union{Symbol, Expr})
    v  = :($s -> $v)
    v1 = :($s -> $v1)
    v2 = :($s -> $v2)
    aexprs = Any[
        Expr(
            :escape,
            Expr(
                :(=),
                Base.Cartesian.inlineanonymous(v, i),
                Expr(
                    :call,
                    :(-),
                    Expr(
                        :call,
                        :(*),
                        Base.Cartesian.inlineanonymous(v1, cycle(i+1)),
                        Base.Cartesian.inlineanonymous(v2, cycle(i+2))
                    ),
                    Expr(
                        :call,
                        :(*),
                        Base.Cartesian.inlineanonymous(v1, cycle(i+2)),
                        Base.Cartesian.inlineanonymous(v2, cycle(i+1))
                    )
                )
            )
        )
        for i = 1:3
    ]
    Expr(:block, aexprs...)
end

#@show @macroexpand  @nexprs(3, u -> $out = @d(v12,u+1)*@d(v32,u+2) - @d(v12,u+2)*@d(v32,u+1))
#@show @macroexpand @nexprs 3 u -> vm_u = @d(v12,u+1)*@d(v32,u+2) - @d(v12,u+2)*@d(v32,u+1)
# f = @macroexpand @cross u out_u a_u b_u
#dump(:(a=11+2))
#dump(f)
# @show f




macro fieldcopy!(dst::Symbol, src::Symbol, components::QuoteNode...)
    ex = quote end
    for comp in components
        push!(ex.args, :(copy!(getproperty($(dst), $(comp)), getproperty($(src), $(comp)))))
    end
    esc(ex)
end

using .XMLRPC: ClientProxy

# macro pymol(ex, obj=nothing, host="localhost", port=9123)
const proxy_defaults = Dict{Symbol,Any}(
    :port=>9123,
    :host=>"http://localhost"
)

export @pymol
macro pymol(ex...)
    body = ex[end]
    kwargs = ex[1:end-1]
    
    obj  = nothing
    port = proxy_defaults[:port]
    host = proxy_defaults[:host]
    
    for kwarg in kwargs
        if Meta.isexpr(kwarg, :(=))
            key,val = kwarg.args
            if isa(key, Symbol)
                if key == :port
                    port = val
                elseif key == :host
                    host = val
                elseif key == :obj
                    obj = val
                end
            else
                throw(ArgumentError("non-symbolic keyword '$key'"))
            end
        else
            throw(ArgumentError("non-keyword argument like option '$kwarg'"))
        end
    end
    return quote
        local val = $(esc(body))
        proxy = getfield($(@__MODULE__), :ClientProxy)($host,$port)
        if val isa State
            proxy.load_coordset(val.coords, $(esc(obj)), 0)
        elseif val isa Pose
            io = IOBuffer()
            write(io, val.graph, val.state)
            proxy.read_pdbstr(String(take!(io)), val.graph.name)
            close(io)
        end
        val
    end
    
end

export @setproperties
macro setproperties(kwargs...)
    target = kwargs[1]
    ex = quote end
    for kwarg in kwargs[2:end]
        if Meta.isexpr(kwarg, :(=))
            key,val = kwarg.args
            if isa(key, Symbol)
                push!(ex.args, :(setproperty!($target, $(Meta.quot(key)), $val)))
            else
                throw(ArgumentError("non-symbolic keyword '$key'"))
            end
        else
            throw(ArgumentError("non-keyword argument like option '$kwarg'"))
        end
    end
    esc(ex)
end

export @isdefined
macro isdefined(var)
    quote
        try
            local _ = $(esc(var))
            true
        catch err
            isa(err, UndefVarError) ? false : rethrow(err)
        end
    end
end