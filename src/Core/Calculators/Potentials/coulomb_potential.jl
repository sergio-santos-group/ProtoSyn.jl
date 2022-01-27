function coulomb_potential(d::T; v::Opt{Tuple{T, T, T}} = nothing, kwargs...) where {T <: AbstractFloat}
    qi = kwargs[:qi]
    qj = kwargs[:qj]

    e = (qi*qj)/d

    v !== nothing && begin
        _f = (qi*qj)/(d * d)
        f1 = - _f .* v
        f2 =   _f .* v
    end

    v !== nothing && return e, f1, f2
    return e
end