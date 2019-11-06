
function equiv(graph::ConnectGraph, atomnums::Vector{Int}, maxdepth::Int=-1)
    
    natoms = length(graph)
    table = zeros(natoms, natoms)
    visited = falses(natoms)
    maxpathlen::Int = 0

    function traverse(pivot::Int, cur::Int, score::Float64=0.0, depth::Int=0)
        
        depth += 1
        if depth > maxpathlen
            maxpathlen = depth
        end

        # increment the score
        #  these magical numbers (0.11 & 0.08) were taken from
        #  http://ambermd.org/antechamber/antechamber.pdf (pag. 4)
        score += (0.11 * depth) + (0.08*atomnums[cur])
        
        # if this is a terminal node, add this path score
        # to the table entry
        if length(graph[cur]) == 1
            table[pivot,cur] += score
        end
        
        if maxdepth > 0 && depth==maxdepth
            return depth
        end

        # mark this node as visited. We will latter mark it as unvisited
        # so that alternative paths involving this node can be found.
        visited[cur] = true
        for i in graph[cur]
            if !visited[i]
                traverse(pivot, i, score, depth)
            end
        end
        visited[cur] = false
        
        
    end

    # fill the equivalence table
    for i = 1:length(graph)
        traverse(i, i)
    end

    # sort along rows
    sort!(table, dims=2)

    # number of terminal atoms (atoms with only a single connection)
    nter = count(l->length(l)==1, values(graph))
    
    
    # to group equivalent atoms, one can do:
    #   d = Dict{UInt,Vector{Int}}()
    #   for (n,r) in enumerate(eachrow(Tst.PS.table))
    #       push!(get!(d, hash(r), Int[]), n)
    #   end

    # return table
    println("maxpathlen = ", maxpathlen)

    table[:, natoms-nter+1:end], maxpathlen

end


# graph = ConnectGraph()

# # test 1
# #  6 <=> 7
# #  3 <=> 5
# # graph[1] = [2]
# # graph[2] = [1,3,4,5]
# # graph[3] = [2,4,6]
# # graph[4] = [2,3,5]
# # graph[5] = [2,4,7]
# # graph[6] = [3]
# # graph[7] = [5]
# # atnums = [1,6,6,6,6,1,1]

# # ANTECHAMBER GRAPH IN FIG 1.
# graph[1]  = [2,5,6,7]
# graph[2]  = [1,3,8,9]
# graph[3]  = [2,4,10,11]
# graph[4]  = [3,12,13,14]
# graph[5]  = [1]
# graph[6]  = [1]
# graph[7]  = [1]
# graph[8]  = [2]
# graph[9]  = [2]
# graph[10] = [3]
# graph[11] = [3]
# graph[12] = [4]
# graph[13] = [4]
# graph[14] = [4]
# atnums = [6,6,6,6,1,1,1,1,1,1,1,1,1,1]

# table = equiv(graph, atnums)



query = ConnectGraph(
    1 => [3],
    2 => [3],
    3 => [1, 2, 4],
    4 => [3],
)

atnums = ["C", "CA", "N", "H"]

function findall(mol::Molecule, query::ConnectGraph, qnames::Vector{String})
    names = [at.name for at in iterbyatom(mol)]
    # names = [at.name
    # for frag in mol.fragments
    #     for at in frag.source.atoms]
            
    table = zeros(Int, length(qnames), length(names))
    
    subgraph = ConnectGraph()
    
    
    for (i,qname) in enumerate(qnames)
        qsize = length(query[i])
        for (j,name) in enumerate(names)
            
            if name==qname && length(mol.bonds[j]) >= qsize
                table[i,j] = 1
                subgraph[j] = [i for i in mol.bonds[j] if names[i] in qnames]
            end
        end
    end

    id2subid = Dict(k => findfirst(table[:,k] .== 1) for k in keys(subgraph))

    name2type = Dict(v=>k for (k,v) in enumerate(unique(names)))
    qtypes = [name2type[q] for q in names]

    reftable = equiv(query, qtypes)
    println(reftable)
    visited = falses(mol.size)
    
    function traverse(cur::Int, path::Vector{Int}, N::Int, fill_counter::Int=0, p2=Int[])
        
        if !haskey(id2subid, cur)
            return
        end

        ptr = id2subid[cur]
        
        push!(path, cur)
        fill_counter += 1

        if fill_counter >= N
            println(path)
        end

        # path[ptr] = cur
        # fill_counter += 1

        # if fill_counter == N
        #     println(" >> ", p2)
        
        #     println(path)
        #     # return
        # end

        visited[cur] = true
        for i in subgraph[cur]
            if !visited[i]# && (haskey(id2subid,i) && !visited[id2subid[i]])
                traverse(i, copy(path), N, fill_counter, copy(p2))
            end
        end
        visited[cur] = false

    end

    qlen = length(qnames)
    for k in keys(subgraph)
        # if id2subid[k] == 1
        #     println(k)
        #     traverse(k, zeros(Int, qlen), qlen)
        # end
        # traverse(k, zeros(Int, qlen), qlen)
        traverse(k, Int[], qlen)
    end
    # uniqnames = unique(names)
    
    # if !issubset(qnames, uniqnames)
    #     return nothing
    # end

    # name2type = Dict(v=>k for (k,v) in enumerate(uniqnames))
    # println(name2type)

    # qtypes = [name2type[name] for name in qnames]

    # query_table, maxdepth = equiv(query, qtypes)
    
    # println(maxdepth)
    # println(query_table)

    # types = [name2type[name] for name in names]
    # query_table, equiv(mol.bonds, types, maxdepth-1)[1], maxdepth
    table, subgraph
end