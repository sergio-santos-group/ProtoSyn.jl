"""
    get_bump_potential_charges(;[c::T = 1.0], [r::T = 1.0]) where {T <: AbstractFloat}

Return a bump potential with charges included, as described in
https://math.stackexchange.com/a/3236066. The bump function is centered around
`c` with a radius `r`.

# See also
[`get_bump_potential`](@ref)

# Examples
```jldoctest
julia> ProtoSyn.Calculators.get_bump_potential_charges(c = 1.0, r = 1.0)
(::ProtoSyn.Calculators.var"#bump_potential_charges#21"{ProtoSyn.Calculators.var"#bump_potential_charges#20#22"{Float64, Float64}}) (generic function with 1 method)
```
"""
function get_bump_potential_charges(A::Union{Type{ProtoSyn.SISD_0}, Type{ProtoSyn.CUDA_2}}; c::T = 1.0, r::T = 1.0) where {T <: AbstractFloat}

    return function bump_potential_charges(d::T; v::Opt{Tuple{T, T, T}} = nothing, qi::T = 0.0, qj::T = 0.0) where {T <: AbstractFloat}

        if -r + c < d < r + c
            a = (d - c) / r
            m = qi * qj

            e = m * exp(-((log(1+a)*log(1+a)) + (log(1-a)*log(1-a))))
        else
            e = 0.0
        end

        v !== nothing && return e, (0.0, 0.0, 0.0), (0.0, 0.0, 0.0)
        return e
    end
end

function get_bump_potential_charges(A::Type{ProtoSyn.SIMD_1}; c::T = 1.0, r::T = 1.0) where {T <: AbstractFloat}

    return function bump_potential_charges(d::T; v::Opt{Vec{4, T}} = nothing, qi::T = 0.0, qj::T = 0.0) where {T <: AbstractFloat}

        if -r + c < d < r + c
            a = (d - c) / r
            m = qi * qj

            e = m * exp(-((log(1+a)*log(1+a)) + (log(1-a)*log(1-a))))
        else
            e = 0.0
        end

        v !== nothing && return e, Vec((0.0, 0.0, 0.0, 0.0)), Vec((0.0, 0.0, 0.0, 0.0))
        return e
    end
end

get_bump_potential_charges(;c::T = 1.0, r::T = 1.0) where {T <: AbstractFloat} = begin
    get_bump_potential_charges(ProtoSyn.acceleration.active, c = c, r = r)
end


"""
    get_bump_potential(;[c::T = 1.0], [r::T = 1.0]) where {T <: AbstractFloat}

Return a bump potential without charges included, as described in
https://math.stackexchange.com/a/3236066. The bump function is centered around
`c` with a radius `r`.

# See also
[`get_bump_potential_charges`](@ref)

# Examples
```jldoctest
julia> ProtoSyn.Calculators.get_bump_potential(c = 1.0, r = 1.0)
(::ProtoSyn.Calculators.var"#bump_potential#25"{ProtoSyn.Calculators.var"#bump_potential#24#26"{Float64, Float64, Float64}}) (generic function with 1 method)
```
"""
function get_bump_potential(;c::T = 1.0, r::T = 1.0, m::T = 1.0) where {T <: AbstractFloat}

    return function bump_potential(d::T; v::Opt{Tuple{T, T, T}} = nothing, kwargs...) where {T <: AbstractFloat}
        if -r + c < d < r + c
            a = (d - c) / r
            e = m * exp(-((log(1+a)*log(1+a)) + (log(1-a)*log(1-a))))
        else
            e = 0.0
        end

        e = (qi*qj)/d

        v !== nothing && begin
            _f = m*exp((-log((d-c)/r+1)^2-log(1-(d-c)/r)^2)*((2*log(1-(d-c)/r))/(r*(1-(d-c)/r))-(2*log((d-c)/r+1))/(r*((d-c)/r+1))))
            f1 = - _f .* v
            f2 =   _f .* v
        end

        v !== nothing && return e, f1, f2
        return e
    end
end