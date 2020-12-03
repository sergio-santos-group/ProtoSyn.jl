module Peptides

# using CUDA
using ..ProtoSyn
using ..ProtoSyn.Builder

export isproline, grammar


# resource directory for this module
const resource_dir = let
    modname = string(nameof(@__MODULE__))
    joinpath(ProtoSyn.resource_dir, modname)
end


include("constants.jl")
include("calculators.jl")
include("Rotamers/Rotamers.jl")


"""
    isproline(r::Residue) -> Bool

Determine if a residue is a proline.
"""
@inline isproline(r::Residue) = uppercase(r.name) == "PRO"



# """
#     peptidejoin(f1::Fragment, f2::Fragment)

# Join fragments `f1` and `f2` through a peptide bond. The bond is made between
# the `C` of the last residue of `f1` and the `N` of the first residue of `f2`.
# `C` will, consequently, be made the parent of `N`; also, the last residue of `f1`
# will be made the parent of the first residue of `f2`.

# # Example
# ```julia-repl
# julia> peptidejoin(frag1, frag2)
# julia> isparent(frag1.graph[end], frag2.graph[1])
# true
# julia> isparent(frag1.graph[end,"C"], frag2.graph[1,"N])
# true
# ```
# """
# peptidejoin(f1::Fragment, f2::Fragment) = begin
#     r1 = f1.graph[end]
#     r2 = f2.graph[1]
#     ProtoSyn.join(r1, "C", r2, "N")
#     # hasparent(r2) && error("r2 is already connected")

#     # atC = r1["C"]
#     # atN = r2["N"]
#     # bond(atC, atN)          # C<->N
#     # setparent!(atN, atC)    # C->N
#     # setparent!(r2, r1)      # r1->r2

#     state = f2.state
    
#     atomstate = state[r2["N"]]
#     @setproperties atomstate b=1.2 θ=deg2rad(120) ϕ=0
    
#     if isproline(r2)
#         atomstate = state[r2["CD"]]
#         @setproperties atomstate θ=deg2rad(122.5) ϕ=0
#     else
#         atomstate = state[r2["H"]]
#         @setproperties atomstate θ=deg2rad(120) ϕ=0
#     end
        
#     atomstate = state[r2["CA"]]
#     @setproperties atomstate θ=deg2rad(120) ϕ=pi
    
#     setoffset!(state, r2["CA"], 0)
#     setoffset!(state, r2["C"],  0)
#     setoffset!(state, r2["O"], pi)
# end


"""
    build_graph!(top::Topology)

Utility function for building the atom and residue graphs of peptides.
The root atom of the atom graph for each segment is taken as the `N` of
the first residue of each segment.
"""
function build_graph!(top::Topology)
    build_tree!(top) do seg
        !isempty(seg) ? get(seg[1].itemsbyname, "N", nothing) : nothing
    end
    top
end


"""
    grammar([T=Float64,] dir::AbstractString=resource_dir)

Build a [`LGrammar`](@ref) for peptides, taking as variables the fragments
in `dir`. The [`peptidejoin`](@ref) function is included as the default operator.
The returned L-grammar is required for building peptides from fragments.

# Examples
```julia-repl
julia> g = Peptides.grammar();
julia> pose = Builder.build(grammar, seq"AAGASTASSE")
...
```
!!! NOTE
    requires update
"""
function grammar(::Type{T}) where {T <: AbstractFloat}
    filename = joinpath(resource_dir, "grammars.yml")
    Builder.fromfile(T, filename, "peptide")
end
grammar() = grammar(Float64)


function build(grammar::LGrammar{T}, derivation, ss::NTuple{3,Number} = SecondaryStructure[:linear]) where {T <: AbstractFloat}

    pose = Builder.build(grammar, derivation)
    Peptides.setss!(pose, ss)
    sync!(pose)
    pose
end


function load(::Type{T}, filename::AbstractString; bonds_by_distance::Bool = false) where {T <: AbstractFloat}

    pose = ProtoSyn.load(T, filename)

    for segment in eachsegment(pose.graph)
        setparent!(segment[1][1], ProtoSyn.origin(pose.graph))

        n_residues = ProtoSyn.count_residues(segment)
        for residue_index in 2:n_residues
            popparent!(segment[residue_index])
            setparent!(segment[residue_index], segment[residue_index - 1])
        end
    end

    if bonds_by_distance
        dm        = ProtoSyn.Calculators.full_distance_matrix(pose)
        threshold = T(0.1)
    end

    atoms   = collect(eachatom(pose.graph))
    n_atoms = length(atoms)
    visited = ProtoSyn.Mask{Atom}(n_atoms)
    for (i, atom_i) in enumerate(atoms)
        visited[atom_i.index] = true
        if bonds_by_distance
            for (j, atom_j) in enumerate(atoms)
                i == j && continue
                atom_j = atoms[j]
                atom_j in atom_i.bonds && continue
                putative_bond = "$(atom_i.symbol)$(atom_j.symbol)"
                !(putative_bond in keys(Peptides.bond_lengths)) && continue
                d = Peptides.bond_lengths[putative_bond]
                d += d * threshold
                dm[i, j] < d && ProtoSyn.bond(atom_i, atom_j)
            end
        end
        for atom_j in atom_i.bonds
            if visited[atom_j.index]
                atom_i.parent = atom_j
            else
                push!(atom_i.children, atom_j)
            end
        end
    end

    reindex(pose.graph)
    ProtoSyn.request_c2i(pose.state)
    sync!(pose)

    pose
end

load(filename::AbstractString; bonds_by_distance::Bool = false) = begin
    Peptides.load(ProtoSyn.Units.defaultFloat, filename, bonds_by_distance = bonds_by_distance)
end

function unbond(pose::Pose, residue_1::Residue, residue_2::Residue)

    isparent(residue_1, residue_2) && return _unbond(pose, residue_1, residue_2)
    isparent(residue_2, residue_1) && return _unbond(pose, residue_2, residue_1)

end

function _unbond(pose::Pose, residue_1::Residue, residue_2::Residue)

    ProtoSyn.unbond(pose, residue_1["C"], residue_2["N"])
    
    # Set correct positioning
    state = pose.state
    _origin = ProtoSyn.origin(pose.graph)
    sync!(pose)

    at_N = state[residue_2["N"]]
    at_N.b = ProtoSyn.distance(at_N, state[_origin])
    at_N.θ = ProtoSyn.angle(at_N, state[_origin], state[_origin.parent])
    v = ProtoSyn.dihedral(at_N, state[_origin], state[_origin.parent], state[_origin.parent.parent])
    at_N.ϕ += v - ProtoSyn.getdihedral(state, residue_2["N"])

    for child in residue_2["N"].children
        at_1 = state[child]
        at_1.θ = ProtoSyn.angle(at_1, at_N, state[_origin])
        v = ProtoSyn.dihedral(at_1, at_N, state[_origin], state[_origin.parent])
        at_1.ϕ += v - ProtoSyn.getdihedral(state, child)

        for grandchild in child.children
            at_2 = state[grandchild]
            v = ProtoSyn.dihedral(at_2, at_1, at_N, state[_origin])
            at_2.ϕ += v - ProtoSyn.getdihedral(state, grandchild)
        end
    end

    ProtoSyn.request_i2c(state, all = true)
end


function uncap!(pose::Pose, residue::Residue)

    @assert length(residue["N"].children) > 2 "Residue $residue doesn't seem to be capped."

    for atom in ProtoSyn.gather(an"H\d+"r(pose), residue)
        pop!(pose, atom)
    end

    # b = 0.9795641888105121
    # θ = 2.095336203538727
    # ϕ = 0.0001555866587524742

    H = Atom!(residue, "H", -3, -1, "H")
    setparent!(H, residue["N"])
    ProtoSyn.bond(residue["N"], H)

    reindex(pose.graph)

    insert!(pose.state, H.index, State(1))
    pose.state[H].b = 0.9795641888105121
    pose.state[H].θ = deg2rad(120)
    pose.state[H].ϕ = deg2rad(180)

    ProtoSyn.request_i2c(pose.state)
end


function rotate_dihedral!(s::State, residue::Residue, dihedral_type::Dihedral.DihedralType, value::T) where {T <: AbstractFloat}

    if dihedral_type == Dihedral.psi
        length(residue.children) == 0 && return s
        residue = residue.children[1]
    end
    ProtoSyn.rotate_dihedral!(s, residue[dihedral_type.atom], value)
    s
end


"""
TO DO
"""
function get_sequence(container::ProtoSyn.AbstractContainer)::String

    sequence = ""
    for residue in eachresidue(container)
        sequence *= three_2_one[residue.name]
    end

    return sequence
end

get_sequence(pose::Pose) = get_sequence(pose.graph)


include("methods.jl")
include("Mutators/Mutators.jl")

end