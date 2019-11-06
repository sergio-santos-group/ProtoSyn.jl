export @ffspec

@doc raw"""

# Examples
ff = ForcefieldSpec()
@ffspec ff :atoms HarmonicBondType ("T1", "T2",  r0,  k)
@ffspec ff :atoms HarmonicBondType ("T1", "T2",  r0,  k) ... ("T1", "T2",  r0,  k)
@ffspec ff :atoms HarmonicBondType begin 
    ("T1", "T2",  r0,  k)
     ...
    ("T1", "T2",  r0,  k)
end

> ff = ForcefieldSpec()
> @macroexpand @ffspec ff.components :bonds HarmonicBondType ("T1", "T2",  1.0,  5000.0)
quote
    container = get!(ff.components, :bonds, Dict())
    container["T1:T2"] = HarmonicBondType("T1:T2", 1.0, 5000.0)
end

> ff = ForcefieldSpec()
> @macroexpand @ffspec ff.components :bonds HarmonicBondType begin
    ("T1", "T2",  1.0,  5000.0),
    ("T3", "T4",  2.0, 10000.0)
end
quote
    container = get!(ff.components, :bonds, Dict())
    container["T1:T2"] = HarmonicBondType("T1:T2", 1.0, 5000.0)
    container["T3:T4"] = HarmonicBondType("T3:T4", 2.0, 10000.0)
end

"""
macro ffspec(container::Expr, comp::QuoteNode, type::Symbol, items...)
    
    if items[1].head == :block
        items = items[1].args[2].args
    end
    
    nf = fieldcount(getfield(Forcefield, type))-2
    ex = quote
        container = get!($(container), $(comp), Dict())
    end

    for item in items
        # get atom type names
        typenames = item.args[1:end-nf-1]
        # get component values (parameters)
        typevalues = item.args[end-nf:end]
        key = genkey(typenames...)
        ex2 = :(container[$key] = $type($key, $(typevalues...)))
        push!(ex.args, ex2)
    end
    esc(ex)
end
