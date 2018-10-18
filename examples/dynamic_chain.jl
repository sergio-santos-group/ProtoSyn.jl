#--------------------------
#      3
#     / \
#    /   \
#   2     4 ..... p[-q-r]-
#   |     |
#   1     5
#
# BONDS:
#   1-2
#   2-3
#   3-4
#   4-5
#   4-p
#
# ANGLES:
#   1-2-3
#   2-3-4
#   3-4-5
#   3-4-p
#   5-4-p
#   4-p-q
#
# DIHEDRALS:
#   1-2-3-4
#   2-3-4-5
#   2-3-4-p
#   3-4-p-q
#   5-4-p-q
#   4-p-q-r

#--------------------------
mutable struct Dihedral
    a1::Int64
    a2::Int64
    a3::Int64
    a4::Int64
end


mutable struct Angle
    a1::Int64
    a2::Int64
    a3::Int64
end
angles = [Angle(1,2,3), Angle(2,3,4),Angle(3,4,5),Angle(3,4,-1),Angle(5,4,-1),Angle(4,-1,-2)]
mutable struct Bond
    a1::Int64
    a2::Int64
end
bonds = [Bond(1,2), Bond(2,3), Bond(3,4), Bond(4,5), Bond(4,-1)]

mutable struct Atom
end

mutable struct ResidueTemplate
    #atoms::Vector{Atom}
    bonds::Vector{Bond}
    angles::Vector{Angle}
    #dihedrals::Vector{Dihedral}
end

mutable struct Residue
    template::ResidueTemplate
    next::Union{Residue,Nothing}
    prev::Union{Residue,Nothing}
end

mutable struct Chain
    
end

println(bonds)
println(angles)
template = ResidueTemplate(bonds, angles)
println(template)