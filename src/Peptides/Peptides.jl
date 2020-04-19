module Peptides

using ..ProtoSyn


# resource directory for this module
const resource_dir = let
    modname = string(nameof(@__MODULE__))
    joinpath(ProtoSyn.resource_dir, modname, "old")
end

include("constants.jl")



@inline isproline(r::Residue) = uppercase(r.name) == "PRO"



function loadresidue(::Type{T}, fname::AbstractString) where {T<:AbstractFloat}
    top, state = read(T, fname, ProtoSyn.PDB)
    
    res, state = pop!(top, state, top[1,1])
    
    #----
    nstate = state[res["N"]]
    nstate.b = T(1.2)
    nstate.θ = T(deg2rad(120))
    nstate.ϕ = T(0)
    # nstate.ϕ = T(pi)

    #----
    if isproline(res)
        xstate = state[res["CD"]]
        angle = T(deg2rad(122.5))
    else
        xstate = state[res["H"]]
        angle = T(deg2rad(120))
    end
    xstate.θ = angle
    xstate.ϕ  = T(0)
    
    #----
    castate = state[res["CA"]]
    castate.θ = angle
    castate.ϕ = T(pi)  

    # setoffset!(state, res["N"],  0)
    setoffset!(state, res["CA"], 0)
    setoffset!(state, res["C"],  0)
    setoffset!(state, res["O"], pi)

    res,state
end

loadresidue(fname::AbstractString) = loadresidue(Float64, fname)



function loaddb(::Type{T}, dir::AbstractString=resource_dir) where {T<:AbstractFloat}
    lib = ResidueDB()
    filenames = filter(f->endswith(f, ".pdb"), readdir(dir))
    for filename in filenames
        r, s = loadresidue(T, joinpath(dir, filename))
        lib[r.name] = (r,s)
    end
    lib
end

loaddb(dir::AbstractString=resource_dir) = loaddb(Float64, dir)


build(::Type{T}, letters::String, db::ResidueDB) where {T} = begin
    build(T, [one_2_three[letter] for letter in letters], db)
end

build(letters::String, db::ResidueDB) = build(Float64, letters, db)

build(seq::Vector{String}, db::ResidueDB)  = build(Float64, seq, db) 

build(::Type{T}, seq::Vector{String}, db::ResidueDB) where {T<:AbstractFloat} = begin
    
    topology = Topology("UNK", 1)   # new topology
    segment = Segment("A", 1)       # add a single chain to the topology
    state = State{T}()              # new state
    push!(topology, segment)
    
    if length(seq) > 0

        append!(segment, state, seq, db, PeptideRxToolbelt)
        
        reindex(topology)

        for atom in eachatom(topology)
            atom.ascendents = ascendents(atom, 4)
        end
        # request internal-to-cartesian conversion
        state.i2c = true
    end
    
    topology, state
end


setoffset!(state::State{T}, at::Atom, default) where T = begin
    # rotates all sibling dihedrals to "at" so that the
    # dihedral angle identified by "at" is equal to "default" 
    if hasparent(at)
        ϕ = state[at].ϕ-default
        for child in at.parent.children
            state[child].ϕ -= ϕ
        end
    end
    state.i2c = true
    state
end



set!(state{T}, residue::Residue, dtype::DihedralTypes.Type, value::T) where {T<:AbstractFloat}= begin
    at = residue[dtype]
    if at !== nothing
        state[at].Δϕ = value
    end
    state
end

setss!(state::State, top::Topology, (ϕ, ψ, ω)::Tuple{Number,Number,Number}) = begin
    t = eltype(state)
    ϕ, ψ, ω = t(ϕ), t(ψ), t(ω)
    for r in eachresidue(top)
        state[r["CA"]].Δϕ = ϕ   # dihedral C-N-CA-C
        state[r["C" ]].Δϕ = ψ   # dihedral N-CA-C-N
        state[r["N" ]].Δϕ = ω   # dihedral CA-C-N-CA
    end
    state.i2c = true
    state
end


#----------------------------------------------------

peptideroot(r::Residue) = get(r, "N")

function peptidejoin(r1::Residue, r2::Residue)
    if hasparent(r2)
        error("r2 is already connected")
    else
        atC = get(r1, "C")
        atN = get(r2, "N")
        push!(atC.bonds, atN)
        push!(atN.bonds, atC)
        # atom graph
        setparent!(atN, atC)
        # residue graph
        setparent!(r2, r1)
    end
end

function peptidesplit(r1::Residue, r2::Residue)
    if !hasparent(r2) || r2.parent !== r1
        error("unable to spli")
    else
        # split residue graph
        popparent!(r2)
        # remove peptide bond
        atC = get(r1, "C")
        atN = get(r2, "N")
        i = findfirst(a->a===atN, atC.bonds)
        j = findfirst(a->a===atC, atN.bonds)
        deleteat!(atC.bonds, i)
        deleteat!(atN.bonds, j)
        popparent!(atN)
    end
end

export PeptideRxToolbelt
const PeptideRxToolbelt = ReactionToolbelt(peptidejoin, peptidesplit, peptideroot)



end