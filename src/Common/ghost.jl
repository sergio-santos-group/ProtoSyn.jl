struct Atom
    index::Int
    name::String
end

struct Residue
    name::String
    atoms::Vector{Atom}
    # bonds::Vector{Int}
    top::Vector{Int}
    size::Int
end

struct Link{T}
    residue1::T
    residue2::T
    # top::
end

struct LinkedResidue
    residue::Residue
    offset::Int
    # linked::Vector{LinkedResidue}
    links::Vector{Link{LinkedResidue}}
end


struct Molecule
    residues::Vector{LinkedResidue}
    links::Vector{Link}
end


aA1 = Atom(1, "A1")
aA2 = Atom(2, "A2")
aA3 = Atom(3, "A3")
aB1 = Atom(1, "B1")
aB2 = Atom(2, "B2")
aB3 = Atom(3, "B3")
aB4 = Atom(4, "B4")


rA = Residue("R1", [aA1, aA2, aA3], [[2],[1,3],[2]], 3)
rB = Residue("R2", [aB1, aB2, aB3, aB4], [[1,3],[1,4],[1,4],[2,3]], 4)

lr1 = LinkedResidue(rA, 0, [])
lr2 = LinkedResidue(rA, 3, [])
lr3 = LinkedResidue(rB, 6, [])

push!(lr1.linked, lr2)
push!(lr2.linked, lr1, lr3)
push!(lr3.linked, lr2)

m = Molecule([lr1, lr2, lr3])



lAA = Link(rA, rA, Dict(:bonds=>[rA.atoms[3],rA.atoms[1]]))


N->C : bondtype1 | _struct, args 
CA->N->C
N->C->CA


function bind(r1::LinkedResidue, r2::LinkedResidue, merge_function::Function)
    ver se r1 já está ligado a r2
    se não:
        criar link(r1,r2)
        adicionar link a r1.links
        adicionar link a r2.links
        merge_function(link)
end

