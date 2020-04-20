module Tst


include("src/ProtoSyn.jl")

using .ProtoSyn
using .ProtoSyn.Peptides

pymol = ProtoSyn.XMLRPC.ServerProxy("http://localhost", 9123)
pymol.delete("all")

function render(t,s,name)
    io = IOBuffer()
    write(io, t, s)
    pymol.read_pdbstr(String(take!(io)), name)
end

lib = ProtoSyn.Peptides.loaddb()
println("done loading DB")
top, state = Peptides.build("AA", lib)
# top, state = Peptides.build("GSPAA", lib)
println("done building")
render(top, state, "init")


if true
    Peptides.setss!(state, top[1], Peptides.SecondaryStructure[:linear])
    # println("done setss!")
    sync!(state, top, true)
    println("done sync!")
    render(top,state,"raw")

    # # render(top, state, "top")

    # # Peptides.setss!(state, top, Peptides.SecondaryStructure[:helix])
    # # sync!(state, top)
    # # render(top, state, "top")



    # # t,s=read("../../covid/obj.pdb", ProtoSyn.PDB)
    # #t,s=read("covid.pdb", ProtoSyn.PDB)
    # #sync!(state, top)
    # #render(t, s, "covid")
    # # println(filter(at->!ProtoSyn.hasparent(at), collect(eachresidue(t))))
    # # println(collect(eachresidue(t)))


    # massage(state::State, top::Topology) = begin
    #     for i=1:state.size
    #         state[i].Δϕ = 0.02*(1-2*rand())
    #     end
    #     ProtoSyn.i2c!(state,top)
    # end

    # #massage(s,t)
    # #render(t, s, "covid")


    # println("start append!")

    # #Peptides.append!(t[1,87], s, ["ALA","GLY"], lib)
    # #reindex(t)
    # #build_tree!(t)
    # #ProtoSyn.i2c!(s,t,true)
    # #render(t, s, "covid2")


    # append!(top[1][end], state, ["SER","GLY","SER"], Tst.lib)
    # insert!(nothing, top[1,1],  state, ["SER","GLY","SER"], Tst.lib)
    # insert!(top[1,2], nothing,  state, ["SER","GLY","SER"], Tst.lib)
    # insert!(top[1,1], top[1,2], state, ["SER","GLY","SER"], Tst.lib)
    
    # Peptides.split(top[1,1], top[1,2])
    # insert!(top[1], state, 1, 2, ["SER","GLY","SER"], lib)

    #append!(top[1][end], state, ["SER","GLY","SER"], lib)
    
    append!(top[1], state, ["SER","GLY","SER"], lib, PeptideRxToolbelt)

    reindex(top)
    for atom in eachatom(top)
        atom.ascendents = ascendents(atom, 4)
    end

    Peptides.setss!(state, top[1], Peptides.SecondaryStructure[:linear])
    ProtoSyn.i2c!(state, top, true)
    render(top,state,"added")

    Peptides.setss!(state, top[1], Peptides.SecondaryStructure[:helix])
    ProtoSyn.i2c!(state, top, true)
    render(top,state,"helix")

end

s = select(top, rname"GLY" & (name"CA" | name"C"))
println(s)

end


# create a new chain in segment (by default is a child of root)
#  - insert!(segment, sequence)
#    <=> append!(segment, sequence)
#
# append a sequence and make it a child o parentID
#  - insert!(segment, parentID, sequence)
#    <=> append!(parentRES, sequence)

# insert a sequence between fromID and toID
#  - insert!(segment, fromID, toID, sequence)
