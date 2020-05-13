module Peptides

using ..ProtoSyn


# resource directory for this module
const resource_dir = let
    modname = string(nameof(@__MODULE__))
    joinpath(ProtoSyn.resource_dir, modname, "old")
end

include("constants.jl")


"""
    isproline(r::Residue) -> Bool

Determine if a residue is a proline.
"""
@inline isproline(r::Residue) = uppercase(r.name) == "PRO"



function reset_ic(pose::Pose{Segment})
    segment = pose.graph
    state = pose.state
    T = eltype(state)

    if length(segment) > 0

        residue = segment[1]
        st = state[residue["N"]]
        st.b = 1.2          # T(1.2)
        st.θ = deg2rad(120) # T(deg2rad(120))
        st.ϕ = 0            # T(0)

        #----
        if isproline(residue)
            st = state[residue["CD"]]
            angle = T(deg2rad(122.5))
        else
            st = state[residue["H"]]
            angle = T(deg2rad(120))
        end
        st.θ = angle
        st.ϕ  = T(0)
        
        #----
        st = state[residue["CA"]]
        st.θ = angle
        st.ϕ = T(pi)  

        setoffset!(state, residue["CA"], 0)
        setoffset!(state, residue["C"],  0)
        setoffset!(state, residue["O"], pi)
    end
    pose
end

"""
    loadfragment[T=Float64,] fname) -> Pose{Fragment,State}

Load a PDB file and convert it into a fragment.
"""
function loadfragment(::Type{T}, fname::AbstractString) where {T<:AbstractFloat}
    # read the full pdb file. The tree has not yet been build
    # because the read function knows nothing about the
    # internal graph structure that we wish to impose
    if endswith(fname, ".pdb")
        pose = read(T, fname, ProtoSyn.PDB)
    else
        error("Unsupported file type: only PDB files are supported")
    end

    # but now we know this is a peptide and can build
    # the graph by passing the appropriate seed finding function
    build_tree!(peptideseed, pose.graph)

    # convert cartesian to internal coordinates (sync)
    ProtoSyn.request_c2i(pose.state; all=true)
    sync!(pose)

    # convert this pose to a fragment and reset internal
    # coordinates
    frag = ProtoSyn.fragment(pose)
    reset_ic(frag)
    
    # rename and return fragment
    frag.graph.name = join(map(r->r.name, eachresidue(frag.graph)), "-")
    frag
end

loadfragment(fname::AbstractString) = loadfragment(Float64, fname)



function loaddb(::Type{T}, dir::AbstractString=resource_dir) where {T<:AbstractFloat}
    lib = ResidueDB()
    foreach(readdir(dir)) do filename
        if endswith(filename, ".pdb")
            frag = loadfragment(T, joinpath(dir, filename))
            lib[frag.graph.name] = frag
        end
    end
    lib
end

loaddb(dir::AbstractString=resource_dir) = loaddb(Float64, dir)


build(::Type{T}, seq::String, db::ResidueDB) where {T<:AbstractFloat} = begin
    
    top = Topology("UNK", 1)
    state = State{T}()
    state.id = top.id
    pose = Pose(top, state)

    if !isempty(seq)
        frag = fragment(T, seq, db)
        frag.graph.name = "A"
        frag.graph.id = 1
        append!(pose, frag, PeptideRxToolbelt)
        
        reindex(top)

        ProtoSyn.request_i2c(state; all=true)
    end

    pose
end

build(letters::String, db::ResidueDB) = build(Float64, letters, db)


"""
build a new fragment
"""
fragment(::Type{T}, seq::String, db::ResidueDB) where {T<:AbstractFloat} = begin
    
    seg = Segment(seq, -1)
    state = State{T}()
    prev = nothing

    for (resid,letter) in enumerate(seq)
        refpose = db[one_2_three[letter]]
        res = copy(refpose.graph)
        st = copy(refpose.state)
        push!(seg, res.items...)
        append!(state, st)
        prev !== nothing && peptidejoin(prev, res)
        prev = res
    end
    reindex(seg)
    seg.id = state.id = ProtoSyn.genid()

    setss!(state, seg, SecondaryStructure[:linear])
    Pose(seg, state)
end


setss!(state::State, seg::Segment, (ϕ, ψ, ω)::NTuple{3,Number}) = begin
    t = eltype(state)
    ϕ, ψ, ω = t(ϕ), t(ψ), t(ω)
    for r in eachresidue(seg)
        setdihedral!(state, r[DihedralTypes.phi], ϕ)
        setdihedral!(state, r[DihedralTypes.psi],  ψ)
        setdihedral!(state, r[DihedralTypes.omega], ω)
        # state[r[@ϕ]].Δϕ = ϕ   # dihedral C-N-CA-C
        # state[r[@ψ]].Δϕ = ψ   # dihedral N-CA-C-N
        # state[r[@ω]].Δϕ = ω   # dihedral CA-C-N-CA
    end
    ProtoSyn.request_i2c(state)
    state
end
setss!(s::State, seg::Segment, (ϕ, ψ)::NTuple{2,Number}) = setss!(s, seg, (ϕ,ψ,pi))


#----------------------------------------------------

peptideroot(r::Residue) = get(r, "N")
peptideroot(s::Segment) = s[1,"N"]

function peptidejoin(r1::Residue, r2::Residue)
    if hasparent(r2)
        error("r2 is already connected")
    else
        atC = get(r1, "C")
        atN = get(r2, "N")
        bond(atC, atN)
        # set graphs
        setparent!(atN, atC)
        setparent!(r2, r1)
    end
end
peptidejoin(s1::Segment, s2::Segment) = peptidejoin(s1[end], s2[1])
peptidejoin(r1::Residue, s2::Segment) = peptidejoin(r1, s2[1])
peptidejoin(s1::Segment, r2::Residue) = peptidejoin(s1[end], r2)

function peptidesplit(r1::Residue, r2::Residue)
    if !hasparent(r2) || r2.parent !== r1
        error("unable to split")
    else
        atC = get(r1, "C")
        atN = get(r2, "N")
        unbond(atC, atN)
        # split graphs
        popparent!(atN)
        popparent!(r2)
    end
end

peptideseed(t::Topology) = [
        s[1,"N"]
        for s in t.items
        if !isempty(s) && haskey(s[1], "N")
    ]


export PeptideRxToolbelt
const PeptideRxToolbelt = ReactionToolbelt(peptidejoin, peptidesplit, peptideroot)



end