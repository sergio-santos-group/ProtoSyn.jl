"""
# TODO
"""
function get_graph_size(atom::Atom; depth::Int = 1, max_depth::Int = 10)
    size  = 1
    depth += 1
    depth > max_depth && return size
    for bond in atom.children
        size += get_graph_size(bond, depth = depth, max_depth = max_depth)
    end

    return size
end


"""
# TODO: Documentation
"""
function sort_children(atom::Atom; rev::Bool = false)

    # 1. Sort by the graph size: small chains first
    atoms = atom.children
    sizes = [get_graph_size(a) for a in atoms]
    order = sortperm(sizes, rev = rev)
    atoms = atoms[order]
    sizes = sizes[order]

    # Case all graph sizes are different
    sizes_set = Set(sizes)
    length(sizes_set) === length(sizes) && return atoms

    # 2. Case multiple chains with the same size: Sort by atom name
    for size in sort(collect(sizes_set))
        sizes_indexes = findall(x -> x === size, sizes)
        v = view(atoms, sizes_indexes)
        sort!(v, by=(a)->a.name, rev = true)
    end

    return atoms
end

function sort_children!(atom::Atom)
    atom.bonds = sort_children(atom)
end


export travel_graph

abstract type SearchAlgorithm end
abstract type BFS <: SearchAlgorithm end
abstract type DFS <: SearchAlgorithm end

"""
# TODO Update Documentation

    travel_graph(start::Atom; [stop::Opt{Atom} = nothing])

Return a `Vector{Atom}` with all atom instances between [`Atom`](@ref) `start`
and `stop`, while following the structure's [Graph](@ref state-types). If no `stop`
[`Atom`](@ref) instance is provided or if it isn't found as a downstream parent
of the `start` [`Atom`](@ref), all instances until no children [`Atom`](@ref)
instances are found are returned (for example, until the end of the current
[Pose](@ref) of [`Segment`](@ref)). Note that the order of the returned
[`Atom`](@ref) instances reflects the organization of the graph followed, and
not the distance/parenthood to the `start` [`Atom`](@ref), and should therefore
be ignored in most cases.

# See also
[`is_contiguous`](@ref) [`hasparent`](@ref) [`setparent!`](@ref)
 
# Examples
```jldoctest
julia> ProtoSyn.travel_graph(pose.graph[1][end][10])
4-element Vector{Atom}:
 Atom{/UNK:1/UNK:1/LEU:21/CD1:334}
 Atom{/UNK:1/UNK:1/LEU:21/HD13:337}
 Atom{/UNK:1/UNK:1/LEU:21/HD12:336}
 Atom{/UNK:1/UNK:1/LEU:21/HD11:335}
```
"""
function travel_graph(start::Atom; stop::Opt{Atom} = nothing, search_algorithm::Type{<: SearchAlgorithm} = ProtoSyn.BFS)
    atoms = Vector{Atom}([start])
    if search_algorithm === ProtoSyn.BFS
        # Note: without `rev = true`, `sort_children` sorts from large to small
        init_bonds = copy(sort_children(start, rev = true))
    else
        init_bonds = copy(sort_children(start))
    end 
    stack = Vector{Atom}(init_bonds)

    while length(stack) > 0
        atom_i = pop!(stack)
        if atom_i == stop
            push!(atoms, atom_i)
            return atoms
        else
            if search_algorithm === ProtoSyn.BFS
                bonds = copy(sort_children(atom_i, rev = true))
            else
                bonds = copy(sort_children(atom_i))
            end 
            stack = vcat(stack, bonds)
        end

        push!(atoms, atom_i)
    end

    return atoms
end