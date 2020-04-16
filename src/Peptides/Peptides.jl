# module Peptides

# using ..ProtoSyn
# using LinearAlgebra: norm, cross

# # resource directory for this module
# const resource_dir = let
#     modname = string(nameof(@__MODULE__))
#     joinpath(ProtoSyn.resource_dir, modname, "old")
# end

# # # pre-defined resources for this module
# # const resources = (
# #     aminoacids = joinpath(resource_dir, "aminoacids.jl"),
# # )

# #include("methods.jl")


# function load(dir::AbstractString=resource_dir)
#     # identify all pdb files in the given directory
#     files = filter(f->endswith(f, ".pdb"), readdir(dir))
    
#     # read molecules and extract residues
#     residues = map(files) do f
#         println(f)
#         mol,state = read(joinpath(dir,f), ProtoSyn.PDB)
#         residue = pop!(mol, mol.residues[1])
#         #residue.metadata = Dict(
#         #    :coords => state.coords
#         #)
#         residue.coords = state.coords
#         residue
#     end

#     # align residues
#     rmat = zeros(3,3)
#     foreach(residues) do residue
#         xyz = residue.coords
#         # xyz = residue.metadata[:coords]
        
#         # identify backbone atoms. The residue will be aligned
#         # such that the N and C atoms lay along the x-axis and the
#         # N, CA, and C atoms in the xy plane
#         idN  = get(residue,  "N").id
#         idCA = get(residue, "CA").id
#         idC  = get(residue,  "C").id

#         # shift coordinates so that atom N is at origin
#         xyz .-= xyz[:, idN]
#         # xyz .-= xyz[idN, :]'

#         # align the N->C vector along x-axis
#         theta = acos(xyz[1,idC] / norm(xyz[:,idC]))
#         axis = cross([1.0, 0.0, 0.0], xyz[:,idC])
#         # theta = acos(xyz[idC,1] / norm(xyz[idC,:]))
#         # axis = cross([1.0, 0.0, 0.0], xyz[idC,:])
#         rotmat!(rmat, axis, -theta)
#         xyz = rmat * xyz
#         # xyz = xyz * rmat'

#         # put the N->CA vector on the xy-plane
#         v = xyz[:,idCA]
#         # v = xyz[idCA,:]
#         theta = sign(v[3])*acos(v[2] / sqrt(v[2]^2 + v[3]^2))
#         rotmat!(rmat, [1.0, 0.0, 0.0], -theta)
        
#         # perform final rotation and update residue coordinates
#         residue.coords = rmat * xyz
#         # residue.metadata[:coords] = rmat * xyz
#         # residue.metadata[:coords] = xyz * rmat'

#     end

#     Dict(r.name=>r for r in residues)
# end

# include("constants.jl")
# include("methods2.jl")

# end



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


function build(::Type{T}, letters::String, db::ResidueDB) where {T}
    build(T, [one_2_three[letter] for letter in letters], db)
end

build(letters::String, db::ResidueDB) = build(Float64, letters, db)

build(seq::Vector{String}, db::ResidueDB)  = build(Float64, seq, db) 
function build(::Type{T}, seq::Vector{String}, db::ResidueDB) where {T<:AbstractFloat}
    
    top = Topology("UNK", 1)    # new topology
    seg = Segment("A", 1)       # add a single chain to the topology
    push!(top, seg)

    state = State{T}()          # new state

    # add residues
    for (i,s) in enumerate(seq)
        println("adding residue $i : $s")
        # copy the residue and its state
        res, st = copy(db[s])

        # add the copy to segment and to the state
        push!(seg, res)
        append!(state, st)
        
        # update the residue ID
        res.id = i
    end

    # renumber topology
    reindex(top)
    println("build: done reindex")

    nres = length(seq)
    if nres > 1
        for i = 2:nres
            # create the petide bond
            atC = get(seg[i-1], "C")
            atN = get(seg[i], "N")
            push!(atC.bonds, atN)
            push!(atN.bonds, atC)
            # push!(atC.node, atN.node)           # atom graph
            # push!(seg[i-1].node, seg[i].node)   # residue graph
            setparent!(atC, atN)           # atom graph
            setparent!(seg[i-1], seg[i])   # residue graph
        end
    end
    if nres > 0
        setparent!(ProtoSyn.origin(top), get(seg[1], "N"))
    end
    for atom in eachatom(top)
        atom.ascendents = ascendents(atom, 4)
    end

    #build_tree!(top)
    #println("build: done build_tree")
    # for res in eachresidue(top)
    #     setoffset!(state, res["N"],  0) # psi
    #     setoffset!(state, res["CA"], 0) # omega
    #     setoffset!(state, res["C"],  0) # phi
    # end

    # request internal-to-cartesian conversion
    #state.i2c = true
    #sync!(state, top, true)
    # i2c!(state, top, true)

    top, state
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

# phi, psi, omega
const SecondaryStructure = Dict{Symbol, Tuple{Number,Number,Number}}()
SecondaryStructure[:helix] = map(deg2rad, (-60.0, -45.0, 180))
SecondaryStructure[:linear] = map(deg2rad, (180.0, 180.0, 180))
SecondaryStructure[:parallel_sheet] = map(deg2rad, (-119.0, 113.0, 180))
SecondaryStructure[:antiparallel_sheet] = map(deg2rad, (-139.0, 135.0, 180))

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



function Base.append!(residue::Residue, state::State, seq::Vector{String}, db::ResidueDB)
    prev = residue

    insertion_point = mapreduce(a->a.index, max, residue.items)
    println("insertion_point: $insertion_point")
    for (i,s) in enumerate(seq)
        # copy the residue and its state
        res, st = copy(db[s])

        # add the copy to segment and to the state
        push!(residue.container, res)
        # append!(state, st)
        insert!(state, insertion_point, st)
        insertion_point += st.size

        # update the residue ID
        res.id = residue.id+i


        #--------------------
        atC = get(prev, "C")
        atN = get(res,  "N")
        push!(atC.bonds, atN)
        push!(atN.bonds, atC)
        setparent!(atC, atN)       # atom graph
        setparent!(prev, res)      # residue graph
        #--------------------

        prev = res
    end
    #reindex(top)

end




end