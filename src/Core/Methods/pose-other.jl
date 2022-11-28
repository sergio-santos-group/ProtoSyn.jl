"""
    recoverfrom!(pose::Pose, backup::Pose)

Recovers the [`Pose`](@ref) `pose` [State](@ref state-types) and [Graph](@ref graph-types) from a
`backup` [`Pose`](@ref) while maintaining any reference to the given
[`Pose`](@ref) `pose`. In essence, when using a `Driver`, simply using `copy!`
will create a new instance, and sometimes this can cause bugs. It's recommended
to employ [`recoverfrom!`](@ref) in such cases.

# Examples
```
julia> ProtoSyn.recoverfrom!(pose, backup)
Pose{Topology}(Topology{/2a3d:31788}, State{Float64}:
 Size: 1140
 i2c: false | c2i: false
 Energy: Dict(:Total => Inf)
)
```
"""
function recoverfrom!(pose::Pose, backup::Pose)
    pose.state = copy(backup.state)
    pose.graph = copy(backup.graph)
    pose.graph.id = backup.graph.id # Set this to maintain ID
    pose.state.id = backup.state.id # Set this to maintain ID
    
    return pose
end


"""
    merge(pose1::Pose, pose2::Pose)

Merge the two given poses, creating a new `Pose` in the process.

# Examples
```
julia> ProtoSyn.merge(pose, pose_mod)
Pose{Topology}(Topology{/merged:32083}, State{Float64}:
 Size: 686
 i2c: false | c2i: false
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
```
julia> ProtoSyn.merge!(pose, pose_mod)
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

Return a symmetry expanded [Pose](@ref pose-types). Create N copies of the given `pose` in
all 3 symmetry axis of a cubic lattice, where `reps` is the number of copies in
each of the dimensions X, Y and Z (N is, therefore, reps[1]*reps[2]*reps[3]).
Length of `reps` must be 3. `unit_cell_dims` sets the distance in each of
dimension to translate the copies, in Angstrom Å. Length of `unit_cell_dims`
must be 3. Copies the given `pose`, returning a new struct.

# See also
`symexp!` `merge`

# Examples
```
julia> ProtoSyn.symexp(pose, [2, 2, 2], [50.0, 50.0, 50.0])
Pose{Topology}(Topology{/UNK:59312}, State{Float64}:
 Size: 9261
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

Return a symmetry expanded [Pose](@ref pose-types). Create N copies of the given `pose` in
all 3 symmetry axis of a cubic lattice, where `reps` is the number of copies in
each of the dimensions X, Y and Z (N is, therefore, reps[1]*reps[2]*reps[3]).
Length of `reps` must be 3. `unit_cell_dims` sets the distance in each of
dimension to translate the copies, in Angstrom Å. Length of `unit_cell_dims`
must be 3. Copies the given `pose`, returning a new struct. Updates/overwrites
the given `pose`.

# See also
`symexp` `merge`

# Examples
```
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