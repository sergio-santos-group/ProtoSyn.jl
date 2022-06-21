"""
    get_graph_size(atom::Atom; depth::Int = 1, max_depth::Int = 10)

Recursivelly search the the [Graph](@ref graph-types) starting from [`Atom`](@ref) `atom`
(inclusive) until no children are identified or the `depth` > `max_depth` (10,
by default).

# Examples
```
julia> ProtoSyn.get_graph_size(pose.graph[1][73]["CA"])
13

julia> ProtoSyn.get_graph_size(pose.graph[1][73]["C"])
3
```
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
    sort_children(atom::Atom; rev::Bool = false)

Sort the given [`Atom`](@ref) `atom` children, by the following criteria:
 1. By [Graph](@ref graph-types) size (follow [Graph](@ref graph-types) by employing the [`get_graph_size`](@ref), small chains first)
 2. By [`Atom`](@ref) name (in case all children chains have the same size; alphabetical order)
By setting `rev` to `true` (`false`, by default), reverses the provided order.

    sort_children!(atom::Atom; rev::Bool = false)

Sort the given [`Atom`](@ref) `atom` children and save the newly sorted childen
in `atom.bonds`.

# Examples
```
julia> ProtoSyn.sort_children!(pose.graph[1][72]["CA"])
3-element Vector{Atom}:
 Atom{/2a3d:31788/A:1/HIS:72/HA:1112}
 Atom{/2a3d:31788/A:1/HIS:72/CB:1113}
 Atom{/2a3d:31788/A:1/HIS:72/C:1124}
```
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

        # Note: Sorting by name is a two step operation.
        # Atoms should be sorted based on the symbol in reverse order (H < C)
        # But equal symbols should be sorted by name in normal order (H1 < H2)

        sort!(v, by=(a) -> a.name, rev = true)
        symbols = collect([a.symbol for a in v])
        length(symbols) === 1 && continue
        symbols_set = collect(Set(symbols))
        for symbol in symbols_set
            symbol_indexes = findall(x -> x === symbol, symbols)
            length(symbol_indexes) === 1 && continue
            v = view(v, symbol_indexes)
            sort!(v, by=(a) -> a.name, rev = false)
        end
    end

    return atoms
end

function sort_children!(atom::Atom; rev::Bool = false)
    atom.bonds = sort_children(atom, rev = rev)
end

# --- TRAVEL GRAPH SECTION -----------------------------------------------------

export travel_graph

abstract type SearchAlgorithm end

struct BFSA <: ProtoSyn.SearchAlgorithm end

function (sa::BFSA)(atom::Atom, stack::Vector{Atom})
    bonds = copy(sort_children(atom))
    @info "$atom -> $bonds"
    return vcat(stack, bonds)
end

"""
    (ProtoSyn.BFS)(atom::Atom, stack::Vector{Atom})

Breath first search algorithm for [`travel_graph`](@ref). Correctly sorts the
given [`Atom`](@ref) `atom` children instances and concatenates with the current
`stack`.

# Examples
```
julia> ProtoSyn.BFS(pose.graph[1][1]["CA"], Vector{Atom}())
3-element Vector{Atom}:
 Atom{/test:36441/A:1/MET:1/HA:6}
 Atom{/test:36441/A:1/MET:1/C:7}
 Atom{/test:36441/A:1/MET:1/CB:8}
```
"""
BFS = BFSA()


struct DFSA <: ProtoSyn.SearchAlgorithm end

function (sa::DFSA)(atom::Atom, stack::Vector{Atom})
    bonds = copy(sort_children(atom, rev = true))
    return vcat(bonds, stack)
end

"""
    (ProtoSyn.DFS)(atom::Atom, stack::Vector{Atom})

Depth first search algorithm for [`travel_graph`](@ref). Correctly sorts the
given [`Atom`](@ref) `atom` children instances and concatenates with the current
`stack`.

# Examples
```
julia> ProtoSyn.DFS(pose.graph[1][1]["CA"], Vector{Atom}())
3-element Vector{Atom}:
 Atom{/test:36441/A:1/MET:1/CB:8}
 Atom{/test:36441/A:1/MET:1/C:7}
 Atom{/test:36441/A:1/MET:1/HA:6}
```
"""
DFS = DFSA()


"""
    travel_graph(start::Atom; [stop::Opt{Atom} = nothing], [search_algorithm::F = ProtoSyn.BFS]) where {F <: SearchAlgorithm})

Return a `Vector{Atom}` with all atom instances between [`Atom`](@ref) `start`
and `stop` (inclusive), while following the structure's
[Graph](@ref state-types). If no `stop` [`Atom`](@ref) instance is provided or
if it isn't found as a downstream parent of the `start` [`Atom`](@ref), all
instances until no children [`Atom`](@ref) instances are found are returned (for
example, until the end of the current [Pose](@ref pose-types) of [`Segment`](@ref)). By
default, uses Breath First Search (BFS) algorithm (all [`Atom`](@ref) instances
at the same "graph-distance" to the `start` [`Atom`](@ref) are consumed before
the next level is considered, order is given by [`sort_children`](@ref)).
Optionally, by setting `search_algorithm` to `ProtoSyn.DFS`, can employ Depth
First Algorithm (DFS) (the largest chain of `atom.children` is recursively
exhausted before consuming the smaller chains, order is given by
[`sort_children`](@ref)).

# See also
[`is_contiguous`](@ref) [`hasparent`](@ref) [`setparent!`](@ref)
 
# Examples
```jldoctest
julia> ProtoSyn.travel_graph(pose.graph[1][5]["N"], stop = pose.graph[1][6]["N"], search_algorithm = ProtoSyn.BFS)
11-element Vector{Atom}:
 Atom{/2a3d:31788/A:1/ALA:5/N:62}
 Atom{/2a3d:31788/A:1/ALA:5/H:63}
 Atom{/2a3d:31788/A:1/ALA:5/CA:64}
 Atom{/2a3d:31788/A:1/ALA:5/HA:65}
 Atom{/2a3d:31788/A:1/ALA:5/CB:66}
 Atom{/2a3d:31788/A:1/ALA:5/C:70}
 Atom{/2a3d:31788/A:1/ALA:5/HB3:69}
 Atom{/2a3d:31788/A:1/ALA:5/HB2:68}
 Atom{/2a3d:31788/A:1/ALA:5/HB1:67}
 Atom{/2a3d:31788/A:1/ALA:5/O:71}
 Atom{/2a3d:31788/A:1/GLU:6/N:72}

julia> ProtoSyn.travel_graph(pose.graph[1][5]["N"], stop = pose.graph[1][6]["N"], search_algorithm = ProtoSyn.DFS)
4-element Vector{Atom}:
 Atom{/2a3d:31788/A:1/ALA:5/N:62}
 Atom{/2a3d:31788/A:1/ALA:5/CA:64}
 Atom{/2a3d:31788/A:1/ALA:5/C:70}
 Atom{/2a3d:31788/A:1/GLU:6/N:72}
```
"""
function travel_graph(start::Atom; stop::Opt{Atom} = nothing, search_algorithm::F = BFS) where {F <: SearchAlgorithm}

    atoms = Vector{Atom}([start])
    stack = search_algorithm(start, Vector{Atom}())

    while length(stack) > 0

        atom_i = stack[1]
        deleteat!(stack, 1)

        if atom_i == stop
            push!(atoms, atom_i)
            return atoms
        else
            stack = search_algorithm(atom_i, stack)
        end

        push!(atoms, atom_i)
    end

    return atoms
end


"""
    travel_graph(start::Atom; [stop::Opt{Atom} = nothing], [search_algorithm::F = ProtoSyn.BFS]) where {F <: SearchAlgorithm})

Return a `Vector{Atom}` with all atom instances between [`Atom`](@ref) `start`
and `stop` (inclusive), while following the structure's
[Graph](@ref state-types). If no `stop` [`Atom`](@ref) instance is provided or
if it isn't found as a downstream parent of the `start` [`Atom`](@ref), all
instances until no children [`Atom`](@ref) instances are found are returned (for
example, until the end of the current [Pose](@ref pose-types) of [`Segment`](@ref)). By
default, uses Breath First Search (BFS) algorithm (all [`Atom`](@ref) instances
at the same "graph-distance" to the `start` [`Atom`](@ref) are consumed before
the next level is considered, order is given by [`sort_children`](@ref)).
Optionally, by setting `search_algorithm` to `ProtoSyn.DFS`, can employ Depth
First Algorithm (DFS) (the largest chain of `atom.children` is recursively
exhausted before consuming the smaller chains, order is given by
[`sort_children`](@ref)).

# See also
[`is_contiguous`](@ref) [`hasparent`](@ref) [`setparent!`](@ref)
 
# Examples
```jldoctest
julia> ProtoSyn.travel_graph(pose.graph[1][5]["N"], stop = pose.graph[1][6]["N"], search_algorithm = ProtoSyn.BFS)
11-element Vector{Atom}:
 Atom{/2a3d:31788/A:1/ALA:5/N:62}
 Atom{/2a3d:31788/A:1/ALA:5/H:63}
 Atom{/2a3d:31788/A:1/ALA:5/CA:64}
 Atom{/2a3d:31788/A:1/ALA:5/HA:65}
 Atom{/2a3d:31788/A:1/ALA:5/CB:66}
 Atom{/2a3d:31788/A:1/ALA:5/C:70}
 Atom{/2a3d:31788/A:1/ALA:5/HB3:69}
 Atom{/2a3d:31788/A:1/ALA:5/HB2:68}
 Atom{/2a3d:31788/A:1/ALA:5/HB1:67}
 Atom{/2a3d:31788/A:1/ALA:5/O:71}
 Atom{/2a3d:31788/A:1/GLU:6/N:72}

julia> ProtoSyn.travel_graph(pose.graph[1][5]["N"], stop = pose.graph[1][6]["N"], search_algorithm = ProtoSyn.DFS)
4-element Vector{Atom}:
 Atom{/2a3d:31788/A:1/ALA:5/N:62}
 Atom{/2a3d:31788/A:1/ALA:5/CA:64}
 Atom{/2a3d:31788/A:1/ALA:5/C:70}
 Atom{/2a3d:31788/A:1/GLU:6/N:72}
```
"""
function travel_bonds(start::Atom, level::Int = 1, previous::Opt{Atom} = nothing)
    if level === 1
        return push!([bond for bond in start.bonds if bond !== previous], start)
    else
        bonds = Vector{Atom}([start])
        for bond in start.bonds
            if bond !== previous
                found_bonds = travel_bonds(bond, level - 1, start)
                bonds = vcat(bonds, found_bonds)
            end
        end
        return bonds
    end
end