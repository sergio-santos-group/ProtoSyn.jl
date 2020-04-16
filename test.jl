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
top, state = Peptides.build("A"^3, lib)
# top, state = Peptides.build("GSPAA", lib)
println("done building")
render(top, state, "init")


if true
    Peptides.setss!(state, top, Peptides.SecondaryStructure[:linear])
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


    append!(top[1][end], state, ["ALA","GLY","SER"], Tst.lib)
    # append!(top[1][end], state, ["ALA"], Tst.lib)
    reindex(top)
    for atom in eachatom(top)
        atom.ascendents = ascendents(atom, 4)
    end
    #     #println(at)
    #     node = at.node
    #     node.ascendents = (
    #         node.item.index,
    #         node.parent.item.index,
    #         node.parent.parent.item.index,
    #         node.parent.parent.parent.item.index
    #     )
    # end
    # # build_tree!(top)

    Peptides.setss!(state, top, Peptides.SecondaryStructure[:linear])
    ProtoSyn.i2c!(state, top, true)
    render(top,state,"added")

    Peptides.setss!(state, top, Peptides.SecondaryStructure[:helix])
    ProtoSyn.i2c!(state, top, true)
    render(top,state,"helix")

end



end
