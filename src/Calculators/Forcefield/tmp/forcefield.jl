
export BondedList, ExclusionList

const BondedList = Dict{Union{Int,Symbol}, Vector{Tuple{Vararg{Int}}}}
const ExclusionList = Dict{Int, Vector{Int}}


function genbonded(mol::Molecule, maxdepth::Int=3)
    if mol.bonds !== nothing
        return genbonded(mol.bonds, maxdepth)
    end
    nothing
end


function genbonded(graph::ConnectGraph, maxdepth::Int=3)
    
    bonded = BondedList()

    visited = falses(length(graph))
    
    function traverse(pivot::Int, path::Vector{Int}, depth::Int=0)
        # do nothing if one is beyond the maximum
        # allowed depth
        if depth > maxdepth
            return
        end
        
        # add this item to the path
        push!(path, pivot)

        # save the path to the adequate container:
        #  by now, the depth correponds to the number
        #  of bonds away from the initial atom 
        if (depth > 0) && (path[1] < path[end])
            push!(get!(bonded, depth, []), Tuple(path))
        end
        
        # if we have reached the maximum depth, then return
        #  and prevent another cycle of traversal that will not
        #  be useful.
        if depth == maxdepth
            return
        end

        # traverse the graph 
        visited[pivot] = true
        for i in graph[pivot]
            if !visited[i]
                traverse(i, copy(path), depth+1)
            end
        end
        visited[pivot] = false
        
    end

    for i = 1:length(graph)
        traverse(i, Int[])
    end

    # if the 1 to 3 keys exist, create named alias
    haskey(bonded, 1) && (bonded[:bonds ] = bonded[1])
    haskey(bonded, 2) && (bonded[:angles] = bonded[2])
    haskey(bonded, 3) && (bonded[:proper] = bonded[3])
    haskey(bonded, maxdepth) && (bonded[:pairs] = bonded[maxdepth])

    return bonded
end





function genexcluded(mol::Molecule, maxdepth::Int=3)
    if mol.bonds !== nothing
        return genexcluded(mol.bonds, maxdepth)
    end
    nothing
end


function genexcluded(graph::ConnectGraph, maxdepth::Int=3)
    
    excluded = ExclusionList()
    visited = falses(length(graph))
    
    function traverse(pivot::Int, path::Vector{Int}, depth::Int=0)
        # do nothing if already beyond the max depth
        (depth > maxdepth) && return
        
        # add this index tp the path only if depth > 0.
        # Otherwise, the first entry of the exclusion list
        # for each atom is the atom itself, and we want the
        # exclusion list to contain only indices greater
        # than the key itsefl
        (depth > 0) && push!(path, pivot)
        
        # prevent the next function call overhead
        (depth == maxdepth) && return

        # mark this atom as visited , traverse the graph
        # and unmark the atom back
        visited[pivot] = true
        for i in graph[pivot]
            !visited[i] && traverse(i,path, depth+1)
        end
        visited[pivot] = false
    end
    

    for k = 1:length(graph)
        # get the list for atom k
        excl = get!(excluded, k, [])
        # build the list and sort it
        traverse(k, excl)
        sort!(excl)
        # mark this atom as visited to prevent
        # indices < k from appearing in the list
        visited[k] = true
    end

    excluded
end
