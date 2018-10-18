mutable struct CallbackObject
    freq::Int64
    f::Function
end

macro callback(freq, f)
    ex = :(CallbackObject($freq, $f))
    println(ex)
    esc(ex)
end

macro faggregator(name, f, args...)
    cachename = Symbol("##", "_faggregator_cache")
    mod = @__MODULE__
    
    cache = isdefined(mod, cachename) ?
             getfield(mod, cachename) :
             Core.eval(mod, :(const $cachename = Dict{Symbol,Expr}()))

    if haskey(cache, name, )
        ex = :(energy += $f($(args...), s, b))
        expr = cache[name]
        pos = length(expr.args[2].args[2].args) - 2
        insert!(expr.args[2].args[2].args, pos, ex)
    else
        cache[name] = quote
            function $name(s::MyState, b::Bool)
                energy = $f($(args...), s, b)
                s.eTotal = energy
                return energy
            end
        end
    end    
    return esc(cache[name])
end

mutable struct MyTop
    x::Int64
end

mutable struct MyState
    eTotal::Int64
    eF::Int64
    eG::Int64
    eH::Int64
end

f(t::MyTop, s::MyState, a...) = (s.eF=100*t.x; s.eF)
g(t::MyTop, s::MyState, a...) = (s.eG= 10*t.x; s.eG)
h(t::MyTop, s::MyState, a...) = (s.eH=  1*t.x; s.eH)

state = MyState(0,0,0,0)
top1 = MyTop(9)
top2 = MyTop(9)
top3 = MyTop(9)

@faggregator myevalf f top1
@faggregator myevalf g top2
@faggregator myevalf h top3

println(myevalf(state, false))
println(state)

#------------
state = MyState(0,0,0,0)
top4 = MyTop(7)
top5 = MyTop(7)

@faggregator myevalg g top4
@faggregator myevalg h top5

println(myevalg(state, false))
println(state)

#-----------------------


a = @callback 10 function(i::Int64, j::Int64)
    i+j
end
dump(a)

# is equivalent to:
function f(i::Int64, j::Int64)
    i+j
end
b = CallbackObject(10, f)
dump(b)

c = CallbackObject(10, (i::Int64, j::Int64)->(i+j))
dump(c)

d = @callback 10 (i::Int64, j::Int64)->(i+j)
dump(d)


function do_work(cb::CallbackObject, n::Int64)
    tmp = cb.f
    for i=1:n
        tmp(1,2)
    end
end

do_work(a, 10)
do_work(b, 10)
do_work(c, 10)
do_work(d, 10)

@time do_work(c, 200000000)
@time do_work(d, 200000000)
@time do_work(a, 200000000)
@time do_work(b, 200000000)
#println(a.f(1,2))


# fb(i,j) = 2*(i+j)
# b = @callback 10 fb
# println(b.f(1,2))

# println(typeof(a))