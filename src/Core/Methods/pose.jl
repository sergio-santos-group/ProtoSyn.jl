"""
    sync!(pose::Pose)
    
Check whether the given `Pose` instance has either i2c or c2i flag set to true
and update the cartesian/internal coordinates accordingly. Return the altered
`Pose` instance.


# Examples
```jldoctest
julia> sync(pose)
```
"""
function sync!(pose::Pose)::Pose
    sync!(pose.state, pose.graph)
    pose
end

function recoverfrom!(pose::Pose, backup::Pose)
    pose.state = copy(backup.state)
    pose.graph = copy(backup.graph)
end