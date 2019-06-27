const Opt = Union{Nothing, T} where T
abstract type AbstractForcefieldComponent end
abstract type AbstractTopologyComponent end


mutable struct Atom
    index::Int64
    name::String
    atomtype::String
    # conects::Vector{Atom}
    residue::String
end

function Atom(input::Vector{Opt{SubString{String}}})
    Atom(
        parse(Int64, input[1]),
        String(input[2]),
        "nothing",
        # Vector{Atom}(),
        String(input[3])
    )
end

function Base.show(io::IO, atom::Atom)
    print(io, "Atom(index=$(atom.index), name=$(atom.name), atomtype=$(atom.atomtype))")
end


# order must match the forcefield (first b, then k)
struct HarmonicBondType <: AbstractForcefieldComponent
    atom_count::Int64
    ordered_atomtypes::Vector{String}
    reverse_ordered_atomtypes::Vector{String}
    b::Float64
    k::Float64
end

function HarmonicBondType(input::Vector{SubString{String}})
    HarmonicBondType(
        2,
        [String(input[1]), String(input[2])], 
        [String(input[2]), String(input[1])],
        parse(Float64, input[3]),
        parse(Float64, input[4])
    )
end


# order must match the forcefield (first b, then k)
struct HarmonicAngleType <: AbstractForcefieldComponent
    atom_count::Int64
    ordered_atomtypes::Vector{String}
    reverse_ordered_atomtypes::Vector{String}
    θ::Float64
    k::Float64
end

function HarmonicAngleType(input::Vector{SubString{String}})
    HarmonicAngleType(
        3,
        [String(input[1]), String(input[2]), String(input[3])], 
        [String(input[3]), String(input[2]), String(input[1])],
        parse(Float64, input[4]),
        parse(Float64, input[5])
    )
end

# -> TODO: ADD DIHEDRAL <-

# order must match the forcefield (first b, then k)
struct HarmonicBond <: AbstractTopologyComponent
    atom1::Atom
    atom2::Atom
    b::Float64
    k::Float64
end


# order must match the forcefield (first b, then k)
struct HarmonicAngle <: AbstractTopologyComponent
    atom1::Atom
    atom2::Atom
    atom3::Atom
    θ::Float64
    k::Float64
end


# name_2_atomtype must be the last item in Forcefield struct
mutable struct Forcefield
    bondtypes::Vector{HarmonicBondType}
    angletypes::Vector{HarmonicAngleType}

    Forcefield() = new(Vector{HarmonicBondType}(), Vector{HarmonicAngleType}())
end


mutable struct Topology
    bonds::Vector{HarmonicBond}
    angles::Vector{HarmonicAngle}

    Topology() = new(Vector{HarmonicBond}(), Vector{HarmonicAngle}())
end