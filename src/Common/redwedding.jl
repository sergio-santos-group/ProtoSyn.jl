const Opt = Union{Nothing, T} where T

struct ResidueMergeRule
    to::Regex
    from::Regex
    allow_intra::Bool
end


mutable struct Atom
    index::Int
    name::String
    degree
end


abstract type AbstractMolecule end

mutable struct Residue{T}
    index::Int
    name::String
    #parent::Opt{AbstractMolecule}
    parent::Opt{T}
    children::Vector{Atom}
    # from::Vector{Residue}
    # to::Vector{Residue}

    #merge_rule::Opt{ResidueMergeRule}
end

# mutable struct Molecule <: AbstractMolecule
mutable struct Molecule
    index::Int
    name::String
    children::Opt{Vector{Residue}}
end


# struct ResiduePack
#     residue::Residue
#     offset::Int
# end

# struct Molecule
#     residues::ResiduePack
# end

# atom() = Atom{Residue}(1,"name",nothing)
# @show a = Atom{Residue}(index=1, name="N")


# 10 ATOMS
atoms = [
    # atoms for residue 1
    Atom(1, "N"),
    Atom(2, "X"),
    Atom(3, "C"),
    # atoms for residue 2
    Atom(4, "N"),
    Atom(5, "C"),
    # atoms for residue 3
    Atom(6, "M"),
    Atom(7, "C"),
    # atoms for residue 4
    Atom(8, "D"),
    Atom(9, "E"),
    Atom(10, "X")
]


# Atom{Residue}(1, "N", nothing)
residues = [
    Residue{Molecule}(1, "VAL", nothing, atoms[1:3]),
    Residue{Molecule}(2, "GLY", nothing, atoms[4:5]),
    Residue{Molecule}(3, "RES", nothing, atoms[6:7]),
    Residue{Molecule}(4, "RES", nothing, atoms[8:10]),
]

# #@show atoms = [Atom(n,string("atom",n)) for n in 1:10]

# # 3 RESIDUES
# @show residues = [Residue(n,string("res",n),nothing,[],[],[],nothing) for n in 1:3]

# @show molecule = Molecule(1,"mol1", residues)

# residues[1].children = atoms[1:3]
# residues[2].children = atoms[4:5]
# residues[3].children = atoms[6:10]

rule  = ResidueMergeRule(r":[A-Z]{3}@(?:(N)):", r":[A-Z]{3}@(?:(C)):", false)
rule2 = ResidueMergeRule(r":[A-Z]{3}@(?:(M)|(D)):", r":[A-Z]{3}@(?:(C)|(E)):", false)


function buildrule(resname::String, from::Vector{String}, to::Vector{String})
    s1 = string("(?:", join([string("(",at,")") for at in from],"|") ,")")
    s2 = string("(?:", join([string("(",at,")") for at in   to],"|") ,")")
    
    ResidueMergeRule(
        Regex(":$(resname)@$(s1):"),
        Regex(":$(resname)@$(s1):"),
        true
    )
end

function buildrule(from::Vector{Pair{String,String}}, to::Vector{Pair{String,String}})

    r1 = join(("($(p.first)@$(p.second))" for p in from), "|")
    r2 = join(("($(p.first)@$(p.second))" for p in   to), "|")

    ResidueMergeRule(Regex(":$r1:"),Regex(":$r2:"),false)
end

buildrule("VAL", ["A1","A2","A3"], ["B1","B2","B3"])

R1@N1 | R1@C3
R2@N2 | R1@C4



function canapply(res1::Residue, res2::Residue, rule::ResidueMergeRule)
    println(res1.index,"<->",res2.index)
    if (res1===res2) && !rule.allow_intra
        return false
    end
        
    names1 = string(":", join(map(at -> res1.name*"@"*at.name, res1.children), ":"), ":")
    names2 = string(":", join(map(at -> res2.name*"@"*at.name, res2.children), ":"), ":")
    
    m1 = match(rule.from, names1)
    m2 = match(rule.to,   names2)
    if !(isa(m1,Nothing) || isa(m2,Nothing))
        println("from> ", names1, m1)
        println("  to> ", names2, m2)
        return true
    end

    m1 = match(rule.from, names2)
    m2 = match(rule.to,   names1)
    if !(isa(m1,Nothing) || isa(m2,Nothing))
        println("from> ", names2, m1)
        println("  to> ", names1, m2)
        return true
    end
    false
end

println(canapply(residues[1], residues[1], rule))
println(canapply(residues[1], residues[2], rule))
println(canapply(residues[1], residues[3], rule))
println(canapply(residues[1], residues[4], rule))

#println(rule)
