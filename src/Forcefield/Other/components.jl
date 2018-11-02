struct SolvPair
    i::Int64
    coef::Float64
end


struct ContactPair
    c1::Int64
    c2::Int64
    prob::Float64
end
Base.show(io::IO, b::ContactPair) = print(io, "Forcefield.Other.ContactPair(c1=$(b.c1), c2=$(b.c2), prob=$(b.prob))")