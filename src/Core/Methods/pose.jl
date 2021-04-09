function recoverfrom!(pose::Pose, backup::Pose)
    pose.state = copy(backup.state)
    pose.graph = copy(backup.graph)
end


"""
    merge(pose1::Pose, pose2::Pose)

Merge the two given poses, creating a new `Pose` in the process.

# Examples
```julia-repl
julia> ProtoSyn.merge(pose1, pose2)
Pose{Topology}(Topology{/merged:10313}, State{Float64}:
 Size: 748
 i2c: false | c2i: true
 Energy: Dict(:Total => Inf)
)
...
```
"""
function merge(pose1::Pose, pose2::Pose)::Pose

    function merge_segments(pose::Pose)
        for segment in pose.graph.items
            s = copy(segment)
            push!(graph, s)
            # Set new parenthood of the first residue in the segment
            ProtoSyn.setparent!(s[1], root.container)
            # Set new parenthood of the first atom in the segment
            ProtoSyn.setparent!(s[1][1], root)
        end
    end

    # Merge graphs
    graph = Topology("merged", -1)
    root = ProtoSyn.root(graph)
    merge_segments(pose1)
    merge_segments(pose2)

    # Merge states
    state = State(pose1.state.size + pose2.state.size)
    state.x[:, 1:pose1.state.size] = pose1.state.x.coords
    state.x[:, (pose1.state.size+1):end] = pose2.state.x.coords

    reindex(graph, set_ascendents = true)
    reindex(state)
    graph.id = state.id = genid()
    sync!(state, graph)
    return Pose(graph, state)
end


"""
    merge!(pose1::Pose, pose2::Pose)

Merge the two given poses, updating/overwritting the given `pose1`.

# Examples
```julia-repl
julia> ProtoSyn.merge!(pose1, pose2)
Pose{Topology}(Topology{/merged:10313}, State{Float64}:
 Size: 748
 i2c: false | c2i: true
 Energy: Dict(:Total => Inf)
)
...
```
"""
function merge!(pose1::Pose, pose2::Pose)

    # Merge graphs
    root = ProtoSyn.root(pose1.graph)
    for segment in pose2.graph.items
        s = copy(segment)
        push!(pose1.graph, s)
        # Set new parenthood of the first residue in the segment
        ProtoSyn.setparent!(s[1], root.container)
        # Set new parenthood of the first atom in the segment
        ProtoSyn.setparent!(s[1][1], root)
    end

    # Merge states (including forces)
    pose1.state.f = hcat(pose1.state.f, pose2.state.f)
    pose1.state.x.coords = hcat(pose1.state.x.coords, pose2.state.x.coords)
    for item in pose2.state.items[4:end]
        item.parent = pose1.state
    end
    pose1.state.items = vcat(pose1.state.items, pose2.state.items[4:end])
    pose1.state.size = length(pose1.state.items) - 3

    reindex(pose1.graph, set_ascendents = true)
    reindex(pose1.state)
    return pose1
end


"""
    symexp(pose::Pose, reps::Vector{Int}, unit_cell_dims::Vector{T}) where {T <: AbstractFloat}

Return a symmetry expanded [Pose](@ref). Create N copies of the given `pose` in
all 3 symmetry axis of a cubic lattice, where `reps` is the number of copies in
each of the dimensions X, Y and Z (N is, therefore, reps[1]*reps[2]*reps[3]).
Length of `reps` must be 3. `unit_cell_dims` sets the distance in each of
dimension to translate the copies, in Angstrom Å. Length of `unit_cell_dims`
must be 3. Copies the given `pose`, returning a new struct.

# See also
`symexp!` `merge`

# Examples
```julia-repl
julia> ProtoSyn.symexp(pose, [2, 2, 2], [50.0, 50.0, 50.0])
Pose{Topology}(Topology{/merged:10313}, State{Float64}:
 Size: 748
 i2c: false | c2i: false
 Energy: Dict(:Total => Inf)
)
...
```
"""
function symexp(pose::Pose, reps::Vector{Int}, unit_cell_dims::Vector{T}) where {T <: AbstractFloat}
    @assert length(reps)==3 "`reps` should be a Vector{Int} with 3 numbers - x, y and z number of repetitions"
    @assert length(unit_cell_dims)==3 "`unit_cell_dims` should be a Vector{AbstractFloat} with 3 numbers - x, y and z distances of the bounding box"

    _pose = copy(pose)
    return symexp!(_pose, reps, unit_cell_dims)
end


"""
    symexp!(pose::Pose, reps::Vector{Int}, unit_cell_dims::Vector{T}) where {T <: AbstractFloat}

Return a symmetry expanded [Pose](@ref). Create N copies of the given `pose` in
all 3 symmetry axis of a cubic lattice, where `reps` is the number of copies in
each of the dimensions X, Y and Z (N is, therefore, reps[1]*reps[2]*reps[3]).
Length of `reps` must be 3. `unit_cell_dims` sets the distance in each of
dimension to translate the copies, in Angstrom Å. Length of `unit_cell_dims`
must be 3. Copies the given `pose`, returning a new struct. Updates/overwrites
the given `pose`.

# See also
`symexp` `merge`

# Examples
```julia-repl
julia> ProtoSyn.symexp!(pose, [2, 2, 2], [50.0, 50.0, 50.0])
Pose{Topology}(Topology{/merged:10313}, State{Float64}:
 Size: 748
 i2c: false | c2i: false
 Energy: Dict(:Total => Inf)
)
...
```
"""
function symexp!(pose::Pose, reps::Vector{Int}, unit_cell_dims::Vector{T}) where {T <: AbstractFloat}
    @assert length(reps)==3 "`reps` should be a Vector{Int} with 3 numbers - x, y and z number of repetitions"
    @assert length(unit_cell_dims)==3 "`unit_cell_dims` should be a Vector{AbstractFloat} with 3 numbers - x, y and z distances of the bounding box"

    x = unit_cell_dims[1]
    y = unit_cell_dims[2]
    z = unit_cell_dims[3]

    i_pose = copy(pose) # pose without the newly added pieces
    for i in 0:reps[1]
        for j in 0:reps[2]
            for k in 0:reps[3]
                _pose = copy(i_pose)
                i==0 && j == 0 && k == 0 && continue
                translation = [i * x, j * y, k * z]
                for i in 1:_pose.state.size
                    _pose.state.x[:, i] = pose.state.x[:, i] .+ translation
                end
                ProtoSyn.request_c2i!(_pose.state)
                sync!(_pose)
                ProtoSyn.merge!(pose, _pose)    
            end # k for
        end # j for
    end # i for
    
    return pose
end # function


export fragment
"""
    fragment(pose::Pose{Topology})
    
Return a [Fragment](@ref) from a given [Pose](@ref) `pose`. The pose must have a
single [`Segment`](@ref).

    fragment(pose::Pose{Topology}, selection::ProtoSyn.AbstractSelection)

Return a [Fragment](@ref) from a list of residues retrieved from the given
`selection` when applied to the provided [Pose](@ref) `pose`. If not yet of
selection type [`Residue`](@ref), the `selection` will be promoted to
[`Residue`](@ref) selection type (with the default `any` aggregating function).
The resulting list of residues must be contiguous (a connected graph of
[`Residue`](@ref) instances parenthoods). These will constitute the unique
[`Segment`](@ref) of the resulting [Fragment](@ref).

!!! ukw "Note:"
    A [Fragment](@ref) is a `Pose{Segment}`, without a root/origin. These are
    usually used as temporary carriers of information, without the ability to be
    directly incorporated in simulations.

# Examples
```jldoctest
julia> frag = fragment(pose)

julia> frag = fragment(pose, rid":1:10")
```
"""
function fragment(pose::Pose{Topology})
    
    length(pose.graph) != 1 && error("only topologies with a single segment can be turned into fragments")
    
    topology = copy(pose.graph)
    segment = topology[1]
    state = splice!(copy(pose.state), 1:count_atoms(segment))
    
    # Detach segment from the old root. This includes removing any parenthood to
    # the origin on any atom or residue.
    detach(segment)

    segment.id = state.id = genid()
    segment.name = topology.name
    segment.container = nothing

    Pose(segment, state)
end


function fragment(pose::Pose{Topology}, selection::ProtoSyn.AbstractSelection)
    # Assumes all residues selected belong to the same Segment

    sele = promote(selection, Residue)
    if !ProtoSyn.is_contiguous(pose, sele)
        error("Tried to fragment a non-contigous selection of residues.")
    end

    @assert length(unique([res.container.id for res in sele(pose, gather = true)])) == 1 "Tried to fragment a contiguous selection of residues belonging to different Segments."

    # Get a copy of the selected residues as a new Segment
    copied_segment = copy(residues[1].container)
    residues       = sele(copied_segment, gather = true)
    segment        = Segment(residues[1].container.name, 1)
    segment.items  = residues
    segment.size   = length(segment.items)

    # Fix the parenthood of the copied residues
    for residue in segment.items
        if !(residue.parent in residues) # parenthood outside this new segment
            residue.parent = nothing
        end
        for child in residue.children
            if !(child in residues) # parenthood to child outside this segment
                indexes = findall(res -> res == child, residue.children)
                deleteat!(residue.children, indexes)
            end
        end

        # Fix the parenthood and bond structures of the copied atoms
        for atom in residue.items
            if !(atom.parent.container in residues)              
                # Remove from bond list
                indexes = findall(at -> at == atom.parent, atom.bonds)
                deleteat!(atom.bonds, indexes)
                
                atom.parent = nothing
            end
            for child in atom.children
                if !(child.container in residues)
                    indexes = findall(at -> at == child, atom.children)
                    deleteat!(atom.children, indexes)

                    # Remove from bond list
                    indexes = findall(at -> at == child, atom.bonds)
                    deleteat!(atom.bonds, indexes)
                end
            end
        end
    end

    # Get a copy of the selected residues' State
    n_atoms = count_atoms(segment)
    state = ProtoSyn.State(n_atoms)
    for (index, atom) in enumerate(eachatom(segment))
        state.items[index + state.index_offset]        = copy(pose.state[atom])
        state.items[index + state.index_offset].parent = state
        state.items[index + state.index_offset].index  = index
        state.x[:, index] = copy(pose.state.x[:, atom.index])
    end

    segment.id = state.id = genid()
    reindex(segment)
    reindex(state)

    return Pose(segment, state)
end