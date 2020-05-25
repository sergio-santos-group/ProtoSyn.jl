
# function loadfragment(::Type{T}, fname::AbstractString, seeder::Function) where {T<:AbstractFloat}
#     # read the full pdb file. The tree has not yet been build
#     # because the read function knows nothing about the
#     # internal graph structure that we wish to impose
#     if endswith(fname, ".pdb")
#         pose = read(T, fname, ProtoSyn.PDB)
#     else
#         error("Unsupported file type: only PDB files are currently supported")
#     end

#     # but now we know this is a peptide and can build
#     # the graph by passing the appropriate seed finding function
#     build_tree!(seeder, pose.graph)

#     # convert cartesian to internal coordinates (sync)
#     request_c2i(pose.state; all=true)
#     sync!(pose)

#     # convert this pose to a fragment and assign a name
#     # based on the input filename
#     frag = fragment(pose)
#     frag.graph.name = uppercase(splitext(basename(fname))[1])

#     frag
# end

# loadfragment(fname::AbstractString, seeder::Function) = loadfragment(Float64, fname, seeder)


# function loaddb(::Type{T}, dir::AbstractString, seeder::Function) where {T<:AbstractFloat}
#     lib = ResidueDB()
#     foreach(readdir(dir)) do filename
#         if endswith(filename, ".pdb")
#             frag = loadfragment(T, joinpath(dir, filename), seeder)
#             lib[frag.graph.name] = frag
#         end
#     end
#     lib
# end

# loaddb(dir::AbstractString, seeder::Function) = loaddb(Float64, dir, seeder)

# function loadresources(::Type{T}, dir::AbstractString) where T
#     resources = []
#     foreach(readdir(dir)) do fname
#         filename = joinpath(dir, fname)
#         if isfile(filename)
#             push!(resources, ProtoSyn.load(T, filename))
#         end 
#     end
#     resources
# end
