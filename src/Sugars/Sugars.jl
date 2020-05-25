module Sugars

using YAML

using ..ProtoSyn
using ..ProtoSyn.Builder


# resource directory for this module
const resource_dir = let
    modname = string(nameof(@__MODULE__))
    joinpath(ProtoSyn.resource_dir, modname)
end


# struct Preset
#     name::AbstractString
#     b::Number
#     θ::Number
#     ϕ::Number
# end

# # function reset_ic(frag::Fragment, presets::Vector{Pair{String,NTuple{3,Number}}})
# function reset_ic(frag::Fragment, presets::Preset...)
#     segment = frag.graph
#     state = frag.state

#     if !isempty(segment)
#         residue = segment[1]
#         for preset in presets
#             atomstate = state[residue[preset.name]]
#             atomstate.θ = deg2rad(preset.θ)
#             atomstate.ϕ = deg2rad(preset.ϕ)
#             if preset.b > 0
#                 atomstate.b = preset.b
#             end
#         end
#     end
#     frag
# end

# const resetmap = Dict(
#     "GLU14" => (Preset("O4", 1.2, 109.5, -62.40), Preset("C4", -1, 120, -180))
#     #"GLU14" => ["C4"=>(120, 120, 180)]
# )
# # const resetmap = Dict(
# #     "GLU14" => (p::Fragment) -> reset_ic(p, "O4"=>(1.2, 120, 0), "C4"=>(-1, 105, 180))
# # )

# # function loadfragment(::Type{T}, fname::AbstractString) where {T<:AbstractFloat}
# #     # read the full pdb file. The tree has not yet been build
# #     # because the read function knows nothing about the
# #     # internal graph structure that we wish to impose
# #     if endswith(fname, ".pdb")
# #         pose = read(T, fname, ProtoSyn.PDB)
# #     else
# #         error("Unsupported file type: only PDB files are supported")
# #     end

# #     # but now we know this is a peptide and can build
# #     # the graph by passing the appropriate seed finding function
# #     build_tree!(sugarseed, pose.graph)

# #     # convert cartesian to internal coordinates (sync)
# #     ProtoSyn.request_c2i(pose.state; all=true)
# #     sync!(pose)

# #     # convert this pose to a fragment and reset internal
# #     # coordinates
# #     frag = ProtoSyn.fragment(pose)
# #     frag.graph.name = splitext(basename(fname))[1]

# #     presets = get(resetmap, frag.graph.name, nothing)
# #     presets!==nothing && reset_ic(frag, presets...)
    
# #     # rename and return fragment
# #     frag
# # end

# # loadfragment(fname::AbstractString) = loadfragment(Float64, fname)

# # function loaddb(::Type{T}, dir::AbstractString=resource_dir) where {T<:AbstractFloat}
# #     lib = ResidueDB()
# #     foreach(readdir(dir)) do filename
# #         if endswith(filename, ".pdb")
# #             frag = loadfragment(T, joinpath(dir, filename))
# #             lib[frag.graph.name] = frag
# #         end
# #     end
# #     lib
# # end

# # loaddb(dir::AbstractString=resource_dir) = loaddb(Float64, dir)

# #------
# function loadfragment(::Type{T}, fname::AbstractString) where {T<:AbstractFloat}
#     frag = ProtoSyn.loadfragment(T, fname, sugarseed)
#     presets = get(resetmap, frag.graph.name, nothing)
#     presets!==nothing && reset_ic(frag, presets...)
#     frag
# end
# loadfragment(fname::AbstractString) = loadfragment(Float64, fname)

# function loaddb(::Type{T}, dir::AbstractString=resource_dir) where {T<:AbstractFloat}
#     lib = ProtoSyn.loaddb(T, dir, sugarseed)
#     for frag in values(lib)
#         presets = get(resetmap, frag.graph.name, nothing)
#         presets!==nothing && reset_ic(frag, presets...)
#     end
#     lib
# end
# loaddb(dir::AbstractString=resource_dir) = loaddb(Float64, dir)


# build(::Type{T}, seq::Vector{String}, db::ResidueDB) where {T<:AbstractFloat} = begin
    
#     top = Topology("UNK", 1)
#     state = State{T}()
#     state.id = top.id
#     pose = Pose(top, state)

#     if !isempty(seq)
#         frag = fragment(T, seq, db)
#         frag.graph.name = "A"
#         frag.graph.id = 1
#         append!(pose, frag, SugarRxToolbelt)
        
#         reindex(top)

#         ProtoSyn.request_i2c(state; all=true)
#     end

#     pose
# end

# build(seq::Vector{String}, db::ResidueDB) = build(Float64, seq, db)

# fragment(::Type{T}, seq::Vector{String}, db::ResidueDB) where {T<:AbstractFloat} = begin
    
#     seg = Segment("A", -1)
#     state = State{T}()
#     prev = nothing

#     for letter in seq
#         refpose = db[letter]
#         res = copy(refpose.graph)
#         st = copy(refpose.state)
#         push!(seg, res.items...)
#         append!(state, st)
#         prev !== nothing && α14(prev, res)
#         prev = res
#     end
#     reindex(seg)
#     seg.id = state.id = ProtoSyn.genid()

#     #setss!(state, seg, SecondaryStructure[:linear])
#     Pose(seg, state)
# end





# sugarseed(t::Topology) = [
#     s[1,"O4"]
#     for s in t.items
#     if !isempty(s) && haskey(s[1], "O4")
# ]
# sugarroot(r::Residue) = get(r, "O4")
# sugarroot(s::Segment) = s[1,"O4"]

# function α14(r1::Residue, r2::Residue)
#     if hasparent(r2)
#         error("r2 is already connected")
#     else
#         atC = get(r1, "C1")
#         atO = get(r2, "O4")
#         bond(atC, atO)
#         # set graphs
#         setparent!(atO, atC)
#         setparent!(r2, r1)
#     end
# end
# α14(s1::Segment, s2::Segment) = α14(s1[end], s2[1])
# α14(r1::Residue, s2::Segment) = α14(r1, s2[1])
# α14(s1::Segment, r2::Residue) = α14(s1[end], r2)

# function sugarsplit(r1::Residue, r2::Residue)
# end

# const SugarRxToolbelt = ReactionToolbelt(α14, sugarsplit, sugarroot)
# end




# function α14(f1::Fragment, f2::Fragment)
#     r1 = f1.graph[end]
#     r2 = f2.graph[1]
#     ProtoSyn.join(r1, "C1", r2, "O4")
    
#     state = f2.state

#     atomstate = state[r2["O4"]]
#     @setproperties atomstate b=1.43 θ=deg2rad(127) ϕ=deg2rad(-90.712)

#     atomstate = state[r2["C4"]]
#     @setproperties atomstate θ=deg2rad(100) ϕ=deg2rad(-103.475)

#     setoffset!(state, r2["H4"], deg2rad(-8.1))


# end

# function build_graph!(top::Topology)
#     build_tree!(top) do seg
#         !isempty(seg) ? get(seg[1].itemsbyname, "O4", nothing) : nothing
#     end
#     top
# end


# function amylose(::Type{T}) where {T<:AbstractFloat}
#     grammar = LGrammar{Char,String}()
#     grammar.defop = α14

#     pose = ProtoSyn.load(T, joinpath(resource_dir, "GLC14.pdb"))
#     if !hasgraph(pose.graph)
#         build_graph!(pose.graph)
#         ProtoSyn.c2i!(pose.state, pose.graph)
#     end
#     grammar['A'] = ProtoSyn.fragment(pose)
#     grammar
# end
# amylose() = amylose(Float64)

function grammar(::Type{T}, key::String) where {T <: AbstractFloat}
    filename = joinpath(resource_dir, "grammars.yml")
    open(filename) do io
        @info "loading grammar from file" filename
        yml = YAML.load(io)
        lgfactory(T, yml[key])
    end
end

grammar(key::String) = grammar(Float64, key)

end
