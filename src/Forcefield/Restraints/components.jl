struct ContactPair
    c1::Int64
    r1::Array{Int64}
    c2::Int64
    r2::Array{Int64}
    prob::Float64
end
Base.show(io::IO, b::ContactPair) = print(io, "Forcefield.Other.ContactPair(c1=$(b.c1), r1=$(b.r1), c2=$(b.c2), r2=$(b.r2), prob=$(b.prob))")