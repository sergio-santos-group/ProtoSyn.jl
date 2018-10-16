macro cbcall(cb, step, args...)
    ex = quote
        if (($cb != nothing) && (getproperty($cb, :freq)>0) && ($step%getproperty($cb, :freq)==0))
            getproperty($cb, :callback)($step, $(args...))
        end
    end
    println(esc(ex))
    return esc(ex)
end