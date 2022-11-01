"""
    sort_atoms_by_graph!(state::State, container::Residue; [start::Opt{Atom} = nothing], [search_algorithm::F = ProtoSyn.BFS]) where {F <: SearchAlgorithm}

Sorts the [`Atom`](@ref) instances in the given [`Residue`](@ref) `container`
to match its graph. By default, employs
[`travel_graph`](@ref ProtoSyn.travel_graph) to get the new sorted list of
[`Atom`](@ref) instances (from the first [`Atom`](@ref) in the
[`Residue`](@ref) `container`, set `start` argument to define a new starting
point). Also updates the [`Atom`](@ref) order in the corresponding and provided
[`State`](@ref) `state`. Expects both the [`State`](@ref) and respective
[Graph](@ref graph-types) to be correctly re-indexed (see
[`reindex`](@ref ProtoSyn.reindex)). By default, uses `search_algorithm` BFS
(breath first search). Note that, after sorting, [`Atom`](@ref) indexes may be
wrong. It's reccommended to [`reindex`](@ref ProtoSyn.reindex) the encompassing
[`Pose`](@ref) after.

    sort_atoms_by_graph!(state::State, container::Union{Topology, Segment}; [start::Opt{Atom} = nothing], [search_algorithm::F = ProtoSyn.BFS]) where {F <: SearchAlgorithm}

Applies [`sort_atoms_by_graph!`](@ref) to all [`Residue`](@ref) instances in the
given `container`. Automatically calls [`reindex`](@ref ProtoSyn.reindex) after
sorting the [`Atom`](@ref) instances.

    sort_atoms_by_graph!(pose::Pose; start::Opt{Atom} = nothing, search_algorithm::F = ProtoSyn.DFS) where {F <: SearchAlgorithm}

Applies [`sort_atoms_by_graph!`](@ref) to all [`Residue`](@ref) instances in the
given [`Pose`](@ref) `pose`. Automatically calls
[`reindex`](@ref ProtoSyn.reindex) after sorting the [`Atom`](@ref) instances.

!!! ukw "Note:"
    When applying [`sort_atoms_by_graph!`](@ref) to a [`Pose`](@ref),
    [`Topology`](@ref) or [`Segment`](@ref), the sorting is still performed at
    the [`Residue`](@ref) level (one [`Residue`](@ref) at a time), therefore the
    chain sorting based on graph size only takes into account intra-residue
    [`Atom`](@ref) instances.

# Examples
```
julia> ProtoSyn.sort_atoms_by_graph!(pose.state, pose.graph[1, 1])
(State{Float64}:
 Size: 1140
 i2c: false | c2i: false
 Energy: Dict(:Total => Inf)
, Residue{/2a3d:3900/A:1/MET:1})

julia> ProtoSyn.sort_atoms_by_graph!(pose.state, pose.graph, search_algorithm = ProtoSyn.Peptides.IUPAC)
(State{Float64}:
 Size: 26
 i2c: false | c2i: false
 Energy: Dict(:Total => Inf)
, Topology{/test:378})
```
"""
function sort_atoms_by_graph!(state::State, container::Residue; start::Opt{Atom} = nothing, search_algorithm::F = ProtoSyn.BFS) where {F <: SearchAlgorithm}

    function find_atom_sortperm(state::State, old_atoms_indexes::Vector{Int}, new_atoms_indexes::Vector{Int})
        _sortperm = Vector{Int}()
        found = 0
        for as in state.items
            if as.index in old_atoms_indexes
                found += 1
                push!(_sortperm, new_atoms_indexes[found])
            else
                push!(_sortperm, as.index)
            end
        end

        @assert found === length(old_atoms_indexes) "From the provided atoms travelled from the graph ($(length(old_atoms_indexes))), only $found were found by the atom.index (in comparison with the atomstate.index). Check if both the graph and state are correctly indexed."

        return _sortperm
    end

    @assert !isa(container, Atom) "Can't sort_atoms_by_graph! on a single Atom instance." 
    old_atoms_indexes = [a.index for a in collect(eachatom(container))]
    
    if start === nothing
        start = ProtoSyn.origin(container)
    end
    if start === nothing
        start = container.items[1]
    end
    @assert start in container "The start atom ($start doesn't belong to the given container ($container)."
    atoms         = ProtoSyn.travel_graph(start, search_algorithm = search_algorithm)
    atoms_indexes = [a.index for a in atoms if a in container.items]

    if length(atoms_indexes) !== length(old_atoms_indexes)
        println([x.name for x in atoms if x in container.items])
        println([x.name for x in collect(eachatom(container))])
    end

    @assert length(atoms_indexes) === length(old_atoms_indexes) "Starting on $start, the inter-container graph traveled only accounts for $(length(atoms_indexes)) of the original $(length(old_atoms_indexes)) atoms. Make sure parenthoods are set and that the defined `start` atom is correct."

    _sortperm = find_atom_sortperm(state, old_atoms_indexes, atoms_indexes)
    @info "Old atom indexes: $old_atoms_indexes"
    @info "New atom indexes: $atoms_indexes"

    # Apply sort perm (ignoring the first 3 entries)
    state.items[4:end] = state.items[4:end][_sortperm[4:end]] # Atom states
    state.x.coords     = state.x.coords[:, _sortperm[4:end]]  # Coords matrix
    container.items    = [a for a in atoms if a in container.items] # Graph

    reindex(state)
    return state, container
end

function sort_atoms_by_graph!(state::State, container::Union{Topology, Segment}; start::Opt{Atom} = nothing, search_algorithm::F = ProtoSyn.DFS, ignore_selection::Opt{AbstractSelection} = nothing) where {F <: SearchAlgorithm}
    
    if ignore_selection === nothing
        ignore_sele = !TrueSelection{Residue}()
    else
        ignore_sele = ProtoSyn.promote(ignore_selection, Residue)
    end
    
    residues = (!ignore_sele)(container, gather = true)
    for residue in residues
        ProtoSyn.sort_atoms_by_graph!(state, residue, start = nothing, search_algorithm = search_algorithm)
    end

    reindex(container)
    return state, container
end

function sort_atoms_by_graph!(pose::Pose; start::Opt{Atom} = nothing, search_algorithm::F = ProtoSyn.DFS, ignore_selection::Opt{AbstractSelection} = nothing) where {F <: SearchAlgorithm}

    if ignore_selection === nothing
        ignore_sele = !TrueSelection{Residue}()
    else
        ignore_sele = ProtoSyn.promote(ignore_selection, Residue)
    end

    return sort_atoms_by_graph!(pose.state, pose.graph, start = start, search_algorithm = search_algorithm, ignore_selection = ignore_sele)
end


include("pose-fragment.jl")
include("pose-modifications.jl")
include("pose-other.jl")
include("pose-diagnose.jl")