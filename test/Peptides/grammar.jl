@testset verbose = true "Peptides | Grammar $(repeat("-", 35))" begin

    for (name, frag) in ProtoSyn.Peptides.grammar.variables
        if isa(frag, Tautomer)
            for tautomer in frag.list
                tname = tautomer.graph[1].name.content
                ts = ProtoSyn.three_2_one[tautomer.graph[1].name.content]
                res = Pose(copy(tautomer)).graph[1, 1]
                graph = ProtoSyn.travel_graph(res[1], sort_bonds = true)
                graph_symbols = [a.symbol for a in graph]
                ProtoSyn.infer_parenthood!(res, overwrite = true, start = res[1])
                new_graph = ProtoSyn.travel_graph(res[1], sort_bonds = true)
                new_graph_symbols = [a.symbol for a in new_graph]
                @test all(graph .=== new_graph)
            end
        else
            res = Pose(copy(frag)).graph[1, 1]
            graph = ProtoSyn.travel_graph(res[1], sort_bonds = true)
            graph_symbols = [a.symbol for a in graph]
            ProtoSyn.infer_parenthood!(res, overwrite = true, start = res[1])
            new_graph = ProtoSyn.travel_graph(res[1], sort_bonds = true)
            new_graph_symbols = [a.symbol for a in new_graph]
            @test all(graph .=== new_graph)
        end
    end
end