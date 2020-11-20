module Peptides

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


"""
    setss!(pose::Pose, (ϕ, ψ, ω)::NTuple{3,Number})

Set the `ϕ`, `ψ` and `ω` backbone angles of all residues in the given `pose`.
This function is usefull for setting the secondary structure of a pose. This
function acts on the internal coordinates and does not update cartesian
coordinates, although a request for conversion is made. It is up to the calling
function/user to explicitly synchornize coordinates via [`sync!`](@ref). In
certain cases it might be useful to not set the origin.children secondary
structure. An example is when a Segment is severed via [`unbond`](@ref), in
which case, updating the origin children will move one of the parts in an big
arm movement.
"""
function setss!(container::Pose, (ϕ, ψ, ω)::NTuple{3,Number}, residues::Vector{Residue})
    state = container.state
    for r in residues
        if r.name == "PRO" # Proline phi is always restricted to -60°
            Builder.setdihedral!(container.state, Dihedral.phi(r), -1.047197551)
        else
            Builder.setdihedral!(container.state, Dihedral.phi(r), ϕ)
        end
        Builder.setdihedral!(container.state, Dihedral.psi(r),  ψ)
        Builder.setdihedral!(container.state, Dihedral.omega(r), ω)
    end
end

setss!(container::Pose, (ϕ, ψ, ω)::NTuple{3, Number}) = begin
    residues::Vector{Residue} = collect(eachresidue(container.graph))
    setss!(container, (ϕ, ψ, ω), residues)
end

setss!(container::Pose, (ϕ, ψ, ω)::NTuple{3,Number}, sele::ProtoSyn.AbstractSelection) = begin
    residues = ProtoSyn.PromoteSelection(sele, Residue, any)(container, gather = true)
    setss!(container, (ϕ, ψ, ω), residues)
end

setss!(container::Pose, (ϕ, ψ, ω)::NTuple{3,Number}, residue::Residue) = begin
    setss!(container, (ϕ, ψ, ω), [residue])
end


function append_residues!(pose::Pose{Topology}, residue::Residue,
    grammar::LGrammar, derivation;
    ss::NTuple{3,Number} = SecondaryStructure[:linear], op = "α")

    Builder.append_residues!(pose, residue, grammar, derivation; op = op)
    residues = residue.container.items[(residue.index + 1):(residue.index + length(derivation))]
    setss!(pose, ss, residues)
end


function insert_residues!(pose::Pose{Topology}, residue::Residue, grammar::LGrammar, derivation;
    ss::NTuple{3, Number} = SecondaryStructure[:linear], op = "α")

    connected_to_origin = residue.parent == ProtoSyn.origin(pose.graph).container
    if connected_to_origin
        state = pose.state
        N = state[residue["N"]] # This is the N atom state
        i = residue.index       # This is the residue index
        (b, θ, ϕ) = (N.b, N.θ, N.ϕ)
    else
        ProtoSyn.unbond(pose, residue.container[residue.index - 1]["C"], residue["N"])
    end
    
    
    # println("Inserting new residues ...")
    Builder.insert_residues!(pose, residue, grammar, derivation; op = op,
        connect_upstream = !connected_to_origin)

    if connected_to_origin
        N = state[residue.container[i]["N"]]
        (N.b, N.θ, N.ϕ) = (b, θ, ϕ)
        ProtoSyn.request_i2c(state; all = true)
    end

    residues = residue.container.items[(residue.index - length(derivation)):(residue.index)]
    setss!(pose, ss, residues)
end


function remove_sidechains!(pose::Pose{Topology}, selection::Opt{AbstractSelection} = nothing)
    _selection = !an"CA$|N$|C$|H$|O$"r
    if selection !== nothing
        _selection = _selection & selection
    end
    sidechain = _selection(pose, gather = true)
    for atom in sidechain
        ProtoSyn.pop_atom!(pose, atom)
    end

    return pose
end


function add_sidechains!(pose::Pose{Topology}, grammar::LGrammar, selection::Opt{AbstractSelection} = nothing)
    if selection !== nothing
        residues = selection(pose, gather = true)
    else
        residues = collect(eachresidue(pose.graph))
    end
    for residue in residues
        derivation = [string(Peptides.three_2_one[residue.name])]
        Peptides.mutate!(pose, residue, grammar, derivation)
    end

    return pose
end


function mutate!(pose::Pose{Topology}, residue::Residue, grammar::LGrammar, derivation)

    @assert length(derivation) == 1 "Derivation must have length = 1."
    
    sidechain = (!an"CA$|N$|C$|H$|O$"r)(residue, gather = true)

    same_aminoacid = string(Peptides.three_2_one[residue.name]) == derivation[1]
    if same_aminoacid && length(sidechain) > 0
        println("No mutation required, residue already has sidechain of the requested type.")
        return pose
    end

    frag = Builder.fragment(grammar, derivation)
    
    # Remove old sidechain
    for atom in sidechain
        ProtoSyn.pop_atom!(pose, atom)
    end
    
    # Insert new sidechain
    frag_sidechain = (!an"CA$|N$|C$|H$|O$"r)(frag, gather = true)
    poseCA = residue["CA"]
    
    # poseCA_index is the LOCAL index (inside residue.items). poseCA.index is
    # the index in the whole pose
    poseCA_index = findfirst(x -> x === poseCA, residue.items)
    
    vals = []
    for (index, atom) in enumerate(frag_sidechain)
        parent_is_CA = false
        if atom.parent.name == "CA"
            parent_is_CA = true
            push!(vals, ProtoSyn.getdihedral(frag.state, atom))
            ProtoSyn.unbond(frag, atom, atom.parent)
        end

        # Insert into the graph
        # Note: insert! already sets the residue.itemsbyname
        insert!(residue, poseCA_index + index, atom)

        # Since now atom is already in the graph, we can bond and add parents
        if parent_is_CA
            ProtoSyn.bond(atom, poseCA)
            setparent!(atom, poseCA)
        end
    end

    _start = frag_sidechain[1].index
    _end   = frag_sidechain[end].index
    insert!(pose.state, poseCA.index + 1, splice!(frag.state, _start:_end))
    
    reindex(pose.graph)
    
    # Fix positions (after reindex - requires correct ascendents)
    pose_sidechain = (!an"CA$|N$|C$|H$|O$"r)(residue, gather = true)
    Δϕ             = pose.state[residue["CA"]].Δϕ
    index          = 1
    for child in residue["CA"].children
        if child in pose_sidechain
            pose.state[child].ϕ = vals[index] - Δϕ
            index += 1
        end
    end
    
    ProtoSyn.request_i2c(pose.state)
    residue.name = one_2_three[derivation[1][1]]

    return pose
end

function pop_residue!(pose::Pose{Topology}, residue::Residue)

    # 1) Since we are calling Peptides, we know that we should unbond this
    # residue and the next (children)
    for child in residue.children
        Peptides.unbond(pose, residue, child) # Order is important
    end

    # 2) Now we can safelly pop the residue, while maintaining the positions of
    # downstream residues
    ProtoSyn.pop_residue!(pose, residue)
end

function build(grammar::LGrammar{T}, derivation, ss::NTuple{3,Number} = SecondaryStructure[:linear]) where {T <: AbstractFloat}

    pose = Builder.build(grammar, derivation)
    setss!(pose, ss)
    sync!(pose)
    pose
end
# build(grammar::LGrammar, derivation, ss::NTuple{3,Number} = SecondaryStructure[:linear]) = build(ProtoSyn.Units.defaultFloat, grammar, derivation, ss)


function load(filename::AbstractString)

    pose = ProtoSyn.load(Float64, filename)

    for segment in eachsegment(pose.graph)
        setparent!(segment[1][1], ProtoSyn.origin(pose.graph))

        n_residues = ProtoSyn.count_residues(segment)
        for residue_index in 2:n_residues
            popparent!(segment[residue_index])
            setparent!(segment[residue_index], segment[residue_index - 1])
        end
    end

    n_atoms = ProtoSyn.count_atoms(pose.graph)
    visited = ProtoSyn.Mask{Atom}(n_atoms)
    for atom_i in eachatom(pose.graph)
        visited[atom_i.index] = true
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


include("Mutators/Mutators.jl")

end