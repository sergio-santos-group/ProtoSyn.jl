"""
    unbond!(pose::Pose, residue_1::Residue, residue_2::Residue; [keep_downstream_position::Bool = false])

Unbond the two provided residues by severing the peptidic bond between the C
atom of `residue_1` and the N atom of `residue_2`. In order to do this, perform
the following steps: unset parent/children, unbond neighbours, remove from
[Graph](@ref graph-types), remove from [State](@ref state-types), update the
containers `itemsbyname` field. If `keep_downstream_position` is set to `true`
(`false` by default), the downstream [`Residue`](@ref) position is maintained
(by calling [`request_c2i!`](@ref ProtoSyn.request_c2i!) and
[`sync!`](@ref ProtoSyn.sync!) methods).

# Examples
```jldoctest
julia> ProtoSyn.Peptides.unbond!(pose, pose.graph[1][2], pose.graph[1][3]; keep_downstream_position = true)
Pose{Topology}(Topology{/UNK:1}, State{Float64}:
 Size: 343
 i2c: false | c2i: false
 Energy: Dict(:Total => Inf)
)
```
"""
function unbond!(pose::Pose, residue_1::Residue, residue_2::Residue; keep_downstream_position::Bool = false)

    if isparent(residue_1, residue_2)
        return ProtoSyn._unbond!(pose, residue_1["C"], residue_2["N"], keep_downstream_position = keep_downstream_position)
    end
    if isparent(residue_2, residue_1)
        return ProtoSyn._unbond!(pose, residue_2["C"], residue_1["N"], keep_downstream_position = keep_downstream_position)
    end
end


"""
    sequence(container::ProtoSyn.AbstractContainer)::String
    sequence(pose::Pose)::String

Return the sequence of aminoacids (in 1 letter mode) of the given container/pose
as a string.

# Examples
```
julia> ProtoSyn.Peptides.sequence(pose)
"SESEAEFKQRLAAIKTRLQAL"
```
"""
function sequence(container::ProtoSyn.AbstractContainer)::String

    sequence = ""
    for residue in eachresidue(container)
        try
            sequence *= three_2_one[residue.name]
        catch KeyError
            sequence *= '?'
        end
    end

    return sequence
end

sequence(pose::Pose) = sequence(pose.graph)