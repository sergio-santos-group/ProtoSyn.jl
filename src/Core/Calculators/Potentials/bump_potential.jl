"""
# TODO: Documentation
"""
function get_bump_potential_charges(;c::T = 1.0, r::T = 1.0) where {T <: AbstractFloat}

    return function bump_potential_charges(d::T; v::Opt{Tuple{T, T, T}} = nothing, qi::T = 0.0, qj::T = 0.0) where {T <: AbstractFloat}

        if -r + c < d < r + c
            a = (d - c) / r
            m = qi * qj

            e = m * exp(-((log(1+a)*log(1+a)) + (log(1-a)*log(1-a))))
            # println("a: $a | m: $m | e: $e")
            if v !== nothing
                _f = m*exp((-log((d-c)/r+1)^2-log(1-(d-c)/r)^2)*((2*log(1-(d-c)/r))/(r*(1-(d-c)/r))-(2*log((d-c)/r+1))/(r*((d-c)/r+1))))
            end
        else
            e = 0.0
            # println("outside zone | e: $e")
            if v !== nothing
                _f = 0.0
            end
        end

        v !== nothing && begin
            f1 = - _f .* v
            f2 =   _f .* v
        end

        v !== nothing && return e, f1, f2
        return e
    end
end


"""
# TODO: Documentation
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