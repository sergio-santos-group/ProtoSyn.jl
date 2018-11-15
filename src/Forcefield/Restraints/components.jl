struct DistanceFBR
    a1::Int64
    a2::Int64
    r1::Float64
    r2::Float64
    r3::Float64
    r4::Float64
    c::Float64
end
Base.show(io::IO, b::DistanceFBR) = print(io, "Forcefield.Restraint.DistanceFBR(a1=$(b.a1), a2=$(b.a2), r1=$(b.r1), r2=$(b.r2), r3=$(b.r3), r4=$(b.r4), c=$(b.c))")


struct DihedralFBR
    a1::Int64
    a2::Int64
    a3::Int64
    a4::Int64
    r1::Float64
    r2::Float64
    r3::Float64
    r4::Float64
    c::Float64
end
Base.show(io::IO, b::DihedralFBR) = print(io, "Forcefield.Restraint.DihedralFBR(a1=$(b.a1), a2=$(b.a2), a3=$(b.a3), a4=$(b.a4), r1=$(b.r1), r2=$(b.r2), r3=$(b.r3), r4=$(b.r4), c=$(b.c))")