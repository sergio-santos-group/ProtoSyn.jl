# module Tst

include("../src/ProtoSyn.jl")

using .ProtoSyn
using .ProtoSyn.Peptides
using .ProtoSyn.Builder

# grammar = Peptides.grammar()
grammar = Peptides.grammar(Float64)
println("Loaded all grammar.")

# @pymol pose = build(grammar, "AAGASTASSE")
pose = build(grammar, seq"AAQG")
println("Done building.")
# println(pose.graph.size)

# Try each one individually
println("\n\nCalculating Saved Mask:")
selection = @atomname "CA"
selection(pose.graph)
@time begin
    selection = @atomname "CA"
    println(selection(pose.graph))
end

println("\n\nUsing Saved Mask:")
@time begin
    println(selection(pose.graph))
end

println("\n\nForcing update:")
@time begin
    println(selection(pose.graph, force_update = true))
end

# selection = @atomsymb "C" | @atomsymb "N"
# println(selection(pose.graph))
# selection = @atomix "H"
# println(selection(pose.graph))
# selection = @teste "ALA"
# println(selection)
# selection2 = @teste "ALA" & selection # DOESN'T WORK
# print(selection2)
# selection = @resname "ALA" | @resname "GLN"
# selection = @resname "ALA" & @resname "GLN"
# selection = (@resname "ALA" | @resname "GLN") & @resname "ALA"

# println(selection(pose.graph))
# exit(1)

# sync!(pose)

# # println(pose.state)

# # ProtoSyn.write(stdout, pose)

# ProtoSyn.write(stdout, pose)

# setss!(pose, SecondaryStructure[:linear])
# @pymol sync!(pose)

# setss!(pose, SecondaryStructure[:helix])
# @pymol sync!(pose)



# end