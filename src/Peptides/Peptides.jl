module Peptides

using YAML

using ..ProtoSyn
using ..ProtoSyn.Builder

export isproline, peptidejoin, grammar, setss!


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
#     loadfragment[T=Float64,] fname) -> Pose{Segment}

# Load a PDB file and convert it into a fragment.
# """
# function peptidesplit(r1::Residue, r2::Residue)
#     if !hasparent(r2) || r2.parent !== r1
#         error("unable to split")
#     else
#         atC = get(r1, "C")
#         atN = get(r2, "N")
#         unbond(atC, atN)
#         # split graphs
#         popparent!(atN)
#         popparent!(r2)
#     end
# end


"""
    peptidejoin(f1::Fragment, f2::Fragment)

Join fragments `f1` and `f2` through a peptide bond. The bond is made between
the `C` of the last residue of `f1` and the `N` of the first residue of `f2`.
`C` will, consequently, be made the parent of `N`; also, the last residue of `f1`
will be made the parent of the first residue of `f2`.

# Example
```julia-repl
julia> peptidejoin(frag1, frag2)
julia> isparent(frag1.graph[end], frag2.graph[1])
true
julia> isparent(frag1.graph[end,"C"], frag2.graph[1,"N])
true
```
"""
peptidejoin(f1::Fragment, f2::Fragment) = begin
    r1 = f1.graph[end]
    r2 = f2.graph[1]
    ProtoSyn.join(r1, "C", r2, "N")
    # hasparent(r2) && error("r2 is already connected")

    # atC = r1["C"]
    # atN = r2["N"]
    # bond(atC, atN)          # C<->N
    # setparent!(atN, atC)    # C->N
    # setparent!(r2, r1)      # r1->r2

    state = f2.state
    
    atomstate = state[r2["N"]]
    @setproperties atomstate b=1.2 θ=deg2rad(120) ϕ=0
    
    if isproline(r2)
        atomstate = state[r2["CD"]]
        @setproperties atomstate θ=deg2rad(122.5) ϕ=0
    else
        atomstate = state[r2["H"]]
        @setproperties atomstate θ=deg2rad(120) ϕ=0
    end
        
    atomstate = state[r2["CA"]]
    @setproperties atomstate θ=deg2rad(120) ϕ=pi
    
    setoffset!(state, r2["CA"], 0)
    setoffset!(state, r2["C"],  0)
    setoffset!(state, r2["O"], pi)
end


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
julia> pose = Builder.build(grammar, "AAGASTASSE")
...
```
"""
function grammar(::Type{T}) where {T <: AbstractFloat}
    filename = joinpath(resource_dir, "grammars.yml")
    open(filename) do io
        @info "loading grammar from file" filename
        yml = YAML.load(io)
        lgfactory(T, yml["peptide"])
    end
end


# function grammar(::Type{T}, dir::AbstractString=resource_dir) where {T<:AbstractFloat}
#     lg = LGrammar{Char,String}()
#     lg.defop = peptidejoin

#     foreach(ProtoSyn.loadresources(T, dir)) do pose
#         if !hasgraph(pose.graph)
#             build_graph!(pose.graph)
#             ProtoSyn.c2i!(pose.state, pose.graph)
#         end
#         frag = ProtoSyn.fragment(pose)
#         println(frag.graph.code)
#         lg[frag.graph.code] = frag
#     end
#     lg
# end
# grammar(dir::AbstractString=resource_dir) = grammar(Float64, dir)


"""
    setss!(pose::Pose, (ϕ, ψ, ω)::NTuple{3,Number})

Set the `ϕ`, `ψ` and `ω` backbone angles of all residues in the given `pose`.
This function is usefull for setting the secondary structure of a pose. This
function acts on the internal coordinates and does not update cartesian
coordinates, although a request for conversion is made. It is up to the calling
function/user to explicitly synchornize coordinates via [`sync!`](@ref). 
"""
setss!(pose::Pose, (ϕ, ψ, ω)::NTuple{3,Number}) = begin
    state = pose.state
    for r in eachresidue(pose.graph)
        setdihedral!(state, r[DihedralTypes.phi], ϕ)
        setdihedral!(state, r[DihedralTypes.psi],  ψ)
        setdihedral!(state, r[DihedralTypes.omega], ω)
    end
    ProtoSyn.request_i2c(state)
end


end