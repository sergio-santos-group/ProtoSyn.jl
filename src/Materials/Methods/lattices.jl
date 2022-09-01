using ProtoSyn

"""
    primitive([unit_cell_dims::Vector{T} = [1.0, 1.0, 1.0]]) where {T <: AbstractFloat}

Creates a primitive lattice `Pose`. If given, the `unit_cell_dims` vector
determines the size of the unit cell. This must be a 1x3 `Vector{T}` of `X`, `Y`
and `Z` lengths of the cubic unit cell (in Angstrom, [1.0, 1.0, 1.0] by default,
using the `ProtoSyn.Units.defaultFloat` as type `T`).

# See also:
[`merge!`](@ref ProtoSyn.merge!) [`symexp!`](@ref ProtoSyn.symexp!)

# Examples
```
julia> pose = ProtoSyn.Materials.primitive()
Pose{Topology}(Topology{/primitive:12909}, State{Float64}:
 Size: 1
 i2c: false | c2i: false
 Energy: Dict(:Total => Inf)
)
```
"""
function primitive(unit_cell_dims::Vector{T} = Vector{ProtoSyn.Units.defaultFloat}([1.0, 1.0, 1.0])) where {T <: AbstractFloat}

    @assert length(unit_cell_dims) == 3 "Unit cell dims must be a 1x3 Vector{T} with X, Y and Z lengths of the cubic unit cell."

    top = Topology("primitive", -1)
    seg = Segment!(top, "unitcell", 1)
    root = ProtoSyn.root(seg)
    res = Residue!(seg, "UNK", 1)
    setparent!(res, root.container)
    atm = Atom!(res, "C", 1, 1, "C")
    setparent!(atm, root)

    state = State(1)
    ProtoSyn.request_c2i!(state)
    top.id = state.id = genid()
    pose = Pose(top, state)
    reindex(pose.graph)
    sync!(pose)

    return pose
end # function


"""
    body_centered([unit_cell_dims::Vector{T} = [1.0, 1.0, 1.0]]) where {T <: AbstractFloat}

Creates a body-centered lattice [`Pose`](@ref). If given, the `unit_cell_dims`
determine the position of the body-centered [`Atom`](@ref) (half of it in all
dimensions). This must be a 1x3 `Vector{T}` of `X`, `Y` and `Z` lengths of the
cubic unit cell (in Angstrom, [1.0, 1.0, 1.0] by default, using the
`ProtoSyn.Units.defaultFloat` as type `T`).

# See also:
[`merge!`](@ref ProtoSyn.merge!) [`symexp!`](@ref ProtoSyn.symexp!)

# Examples
```
julia> ProtoSyn.Materials.body_centered()
Pose{Topology}(Topology{/primitive:51728}, State{Float64}:
    Size: 2
    i2c: false | c2i: false
    Energy: Dict(:Total => Inf)
)
```
"""
function body_centered(unit_cell_dims::Vector{T} = Vector{ProtoSyn.Units.defaultFloat}([1.0, 1.0, 1.0])) where {T <: AbstractFloat}

    @assert length(unit_cell_dims) == 3 "Unit cell dims must be a 1x3 Vector{T} with X, Y and Z lengths of the cubic unit cell."

    top = Topology("body_centered", -1)
    seg = Segment!(top, "unitcell", 1)
    root = ProtoSyn.root(seg)
    res = Residue!(seg, "UNK", 1)
    setparent!(res, root.container)
    atm1 = Atom!(res, "C1", 1, 1, "C")
    setparent!(atm1, root)
    atm2 = Atom!(res, "C2", 1, 1, "C")
    setparent!(atm2, atm1)

    state = State(2)
    x2 = unit_cell_dims[1]/2
    y2 = unit_cell_dims[2]/2
    z2 = unit_cell_dims[3]/2
    state.items[5].t = Vector{T}([x2, y2, z2])
    ProtoSyn.request_c2i!(state)
    top.id = state.id = genid()
    pose = Pose(top, state)
    reindex(pose.graph)
    sync!(pose)

    return pose
end # function


"""
    face_centered([unit_cell_dims::Vector{T} = [1.0, 1.0, 1.0]]) where {T <: AbstractFloat}

Creates a face-centered lattice [`Pose`](@ref). If given, the `unit_cell_dims`
determine the position of the face-centered [`Atom`](@ref) instnaces (half of it
in all dimensions). This must be a 1x3 `Vector{T}` of `X`, `Y` and `Z` lengths
of the cubic unit cell (in Angstrom, [1.0, 1.0, 1.0] by default, using the
`ProtoSyn.Units.defaultFloat` as type `T`).

# See also:
[`merge!`](@ref ProtoSyn.merge!) [`symexp!`](@ref ProtoSyn.symexp!)

# Examples
```
julia> ProtoSyn.Materials.face_centered()
Pose{Topology}(Topology{/face_centered:11920}, State{Float64}:
 Size: 4
 i2c: false | c2i: false
 Energy: Dict(:Total => Inf)
)
```
"""
function face_centered(unit_cell_dims::Vector{T} = Vector{ProtoSyn.Units.defaultFloat}([1.0, 1.0, 1.0])) where {T <: AbstractFloat}

    @assert length(unit_cell_dims) == 3 "Unit cell dims must be a 1x3 Vector{T} with X, Y and Z lengths of the cubic unit cell."

    top = Topology("face_centered", -1)
    seg = Segment!(top, "unitcell", 1)
    root = ProtoSyn.root(seg)
    res = Residue!(seg, "UNK", 1)
    setparent!(res, root.container)
    atm1 = Atom!(res, "C1", 1, 1, "C")
    setparent!(atm1, root)
    atm2 = Atom!(res, "C2", 1, 1, "C")
    setparent!(atm2, atm1)
    atm3 = Atom!(res, "C3", 1, 1, "C")
    setparent!(atm3, atm1)
    atm4 = Atom!(res, "C4", 1, 1, "C")
    setparent!(atm4, atm1)

    state = State(4)
    x2 = unit_cell_dims[1]/2
    y2 = unit_cell_dims[2]/2
    z2 = unit_cell_dims[3]/2
    state.items[5].t = Vector{T}([0, y2, z2]) # atm1
    state.items[6].t = Vector{T}([x2, 0, z2]) # atm2
    state.items[7].t = Vector{T}([x2, y2, 0]) # atm3
    ProtoSyn.request_c2i!(state)
    top.id = state.id = genid()
    pose = Pose(top, state)
    reindex(pose.graph)
    sync!(pose)

    return pose
end # function