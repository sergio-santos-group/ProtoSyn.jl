module Peptides

using ..ProtoSyn
using ..ProtoSyn.Builder

export isproline, grammar, setss!


# resource directory for this module
const resource_dir = let
    modname = string(nameof(@__MODULE__))
    joinpath(ProtoSyn.resource_dir, modname)
end


include("constants.jl")


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
arm movement. This can be controled with the flag 'include_origin_children'
(true by default).
"""
function setss!(container::Pose, (ϕ, ψ, ω)::NTuple{3,Number}, residues::Vector{Residue}; include_origin_children = true)
    state = container.state
    for r in residues
        if !include_origin_children & any(x -> x in ProtoSyn.origin(r).children, r.items)
            continue
        end
        setdihedral!(container, r, Dihedral.phi, ϕ)
        setdihedral!(container, r, Dihedral.psi,  ψ)
        setdihedral!(container, r, Dihedral.omega, ω)
    end
    ProtoSyn.request_i2c(state)
end

setss!(container::Pose, (ϕ, ψ, ω)::NTuple{3, Number}; include_origin_children = true) = begin
    residues::Vector{Residue} = collect(eachresidue(container.graph))
    setss!(container, (ϕ, ψ, ω), residues, include_origin_children = include_origin_children)
end

setss!(container::Pose, (ϕ, ψ, ω)::NTuple{3,Number}, sele::ProtoSyn.AbstractSelection; include_origin_children = true) = begin
    residues = ProtoSyn.PromoteSelection(sele, Residue, any)(container, gather = true)
    setss!(container, (ϕ, ψ, ω), residues, include_origin_children = include_origin_children)
end

setss!(container::Pose, (ϕ, ψ, ω)::NTuple{3,Number}, residue::Residue; include_origin_children = true) = begin
    setss!(container, (ϕ, ψ, ω), [residue], include_origin_children = include_origin_children)
end


function append_residues!(pose::Pose{Topology}, residue::Residue,
    grammar::LGrammar, derivation;
    ss::NTuple{3,Number} = SecondaryStructure[:linear], op = "α")

    Builder.append_residues!(pose, residue, grammar, derivation; op = op)
    residues = residue.container.items[end-length(derivation)+1:end]
    setss!(pose, ss, residues)
end


function insert_residues!(pose::Pose{Topology}, residue::Residue, grammar::LGrammar, derivation;
    ss::NTuple{3, Number} = SecondaryStructure[:linear], op = "α")

    connected_to_origin = residue.parent == ProtoSyn.origin(pose.graph).container
    connect_upstream = false
    if connected_to_origin
        state = pose.state
        N = state[residue["N"]] # This is the N atom state
        i = residue.index       # This is the residue index
        (b, θ, ϕ) = (N.b, N.θ, N.ϕ)
    else
        unbond(pose, residue.container[residue.index - 1]["C"], residue["N"])
        connect_upstream = true
    end
    
    
    # println("Inserting new residues ...")
    Builder.insert_residues!(pose, residue, grammar, derivation; op = op, connect_upstream = connect_upstream)

    if connected_to_origin
        N = state[residue.container[i]["N"]]
        (N.b, N.θ, N.ϕ) = (b, θ, ϕ)
        ProtoSyn.request_i2c(state; all = true)
    end

    residues = residue.container.items[(residue.index - length(derivation)):(residue.index)]
    setss!(pose, ss, residues)
end


function mutate!(pose::Pose{Topology}, residue::Residue, grammar::LGrammar, derivation; op = "α")

    @assert length(derivation) == 1 "Derivation must have length = 1."
    
    frag = Builder.fragment(grammar, derivation)
    measure_frag = Builder.fragment(grammar, vcat(derivation, derivation))

    # Measure angle between Psi and CB_dihedral
    frag_psi = measure_frag.state[measure_frag.graph[2][Dihedral.psi.atom]].ϕ
    measure_frag_sidechain = (!an"CA$|N$|C$|H$|O$"r)(measure_frag, gather = true)
    vals = []
    for child in measure_frag.graph[1]["CA"].children
        if child in measure_frag_sidechain
            push!(vals, frag_psi - measure_frag.state[child].ϕ)
        end
    end
    # frag_CBϕ = measure_frag.state[measure_frag.graph[1]["CB"]].ϕ
    # println("FRAG PSI: $(rad2deg(frag_psi))")
    for item in vals
        println("ITEM: $(rad2deg(item))")
    end
    
    # Remove old sidechain
    sidechain = (!an"CA$|N$|C$|H$|O$"r)(residue, gather = true)
    for atom in sidechain
        pop!(pose, atom)
    end

    # Insert new sidechain
    frag_sidechain = (!an"CA$|N$|C$|H$|O$"r)(frag, gather = true)
    poseCA = residue["CA"]
    poseCA_index = findfirst(x -> x === poseCA, residue.items)
        
    fragCA_children_names = []
    for (index, atom) in enumerate(frag_sidechain)
        if atom.parent.name == "CA"
            push!(fragCA_children_names, atom.name)
            ProtoSyn.unbond(pose, atom, atom.parent)
        end
        insert!(residue, poseCA_index + index, atom)
        residue.itemsbyname[atom.name] = atom
    end

    _start = frag_sidechain[1].index
    _end   = frag_sidechain[end].index
    insert!(pose.state, poseCA.index + 1, splice!(frag.state, _start:_end))
    
    fragCA = frag.graph[1]["CA"]
    index  = 1
    for name in fragCA_children_names

        ProtoSyn.bond(residue[name], poseCA)

        setparent!(residue[name], poseCA)

    end
    
    reindex(pose.graph)
    
    # Fix positions (after reindex - requires correct ascendents ?)
    # Measure angle between Psi and CB_dihedral
    pose_psi = pose.state[residue.container[residue.index + 1][Dihedral.psi.atom]].ϕ
    pose_sidechain = (!an"CA$|N$|C$|H$|O$"r)(residue, gather = true)
    index = 1
    println(residue["CA"].children)
    for child in residue["CA"].children
        if child in pose_sidechain
            pose.state[child].ϕ = (pose_psi + vals[index]) % 360
            index += 1
        end
    end
    
    ProtoSyn.request_i2c(pose.state)
    residue.name = one_2_three[derivation[1][1]]

    return pose
end


function setdihedral!(pose::Pose{Topology}, residue::Residue, dihedral_type::Dihedral.DihedralType, value::T) where {T <: AbstractFloat}
    ProtoSyn.setdihedral!(pose.state, residue[dihedral_type.atom], value)
    pose
end


function getdihedral(pose::Pose{Topology}, residue::Residue, dihedral_type::Dihedral.DihedralType)
    # if dihedral_type.atom == Dihedral.psi.atom
    #     println("CALLED")
    #     length(residue.container.items) <= residue.index && return
    #     residue = residue.container[residue.index + 1]
    # end
    println("Looking for dihedral $(dihedral_type.name) in atom $(dihedral_type.atom)-$(residue[dihedral_type.atom]) (ASCEDENTS: $(residue[dihedral_type.atom].ascendents))")
    r = ProtoSyn.getdihedral!(pose.state, residue[dihedral_type.atom])
    println("  $(rad2deg(r))")
    r
end


function build(::Type{T}, grammar::LGrammar, derivation,
    ss::NTuple{3,Number} = SecondaryStructure[:linear]) where {T<:AbstractFloat}

    pose = Builder.build(T, grammar, derivation)
    setss!(pose, ss, include_origin_children = false)
    pose
end
build(grammar::LGrammar, derivation, ss::NTuple{3,Number} = SecondaryStructure[:linear]) = build(Float64, grammar, derivation, ss)


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

end