"""
    module Lattices

Makes available methods for the initialization of certain lattice poses (such as
`primitive` and `body_centered`).

# See also:
`merge!` `symexp!`
"""
module Lattices

    using ProtoSyn

    """
        primitive([unit_cell_dims::Vector{T} = [1.0, 1.0, 1.0]]) where {T <: AbstractFloat}

    Creates a primitive lattice `Pose`.

    # Examples
    ```jldoctest
    julia> pose = ProtoSyn.Materials.Lattices.primitive()
    Pose{Topology}(Topology{/primitive:25216}, State{Float64}:
     Size: 1
     i2c: false | c2i: false
     Energy: Dict(:Total => Inf)
    )
    ```

    # See also:
    `merge!` `symexp!`
    """
    function primitive(unit_cell_dims::Vector{T} = Vector{ProtoSyn.Units.defaultFloat}([1.0, 1.0, 1.0])) where {T <: AbstractFloat}

        @assert length(unit_cell_dims) == 3 "Unit cell dims must be a 1x3 Vector{T} with X, Y and Z lengths of the cubic unit cell."

        top = Topology("primitive", -1)
        seg = Segment!(top, "unitcell", 1)
        root = ProtoSyn.origin(seg)
        res = Residue!(seg, "UNK", 1)
        setparent!(res, root.container)
        atm = Atom!(res, "C", 1, 1, "C")
        setparent!(atm, root)

        state = State(1)
        ProtoSyn.request_c2i(state)
        top.id = state.id = genid()
        pose = Pose(top, state)
        reindex(pose.graph)
        sync!(pose)

        return pose
    end # function


    """
        body_centered([unit_cell_dims::Vector{T} = [1.0, 1.0, 1.0]]) where {T <: AbstractFloat}

    Creates a body-centered lattice `Pose`. If given, the `unit_cell_dims`
    determine the position of the body-centered atom (half of it in all
    dimensions). This must be a 1x3 Vector{T} of X, Y and Z lengths of the cubic
    unit cell (in Angstrom).

    # Examples
    ```jldoctest
    julia> ProtoSyn.Materials.Lattices.body_centered()
    Pose{Topology}(Topology{/primitive:51728}, State{Float64}:
     Size: 2
     i2c: false | c2i: false
     Energy: Dict(:Total => Inf)
    )
    ```

    # See also:
    `merge!` `symexp!`
    """
    function body_centered(unit_cell_dims::Vector{T} = Vector{ProtoSyn.Units.defaultFloat}([1.0, 1.0, 1.0])) where {T <: AbstractFloat}

        @assert length(unit_cell_dims) == 3 "Unit cell dims must be a 1x3 Vector{T} with X, Y and Z lengths of the cubic unit cell."

        top = Topology("primitive", -1)
        seg = Segment!(top, "unitcell", 1)
        root = ProtoSyn.origin(seg)
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
        ProtoSyn.request_c2i(state)
        top.id = state.id = genid()
        pose = Pose(top, state)
        reindex(pose.graph)
        sync!(pose)

        return pose
    end # function

end # Lattices module