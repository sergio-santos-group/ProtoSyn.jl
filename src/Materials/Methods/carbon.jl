function new_layer(x::Int, y::Int, r::T = 1.4) where {T <: AbstractFloat}

    state    = State()
    topol    = Topology("CRV", -1)
    topol.id = state.id = genid()
    pose     = Pose(topol, state)
    root     = ProtoSyn.root(topol)
    Segment!(topol, "CRV", -1)
    Residue!(topol[1], "CRV", -1)
    ProtoSyn.setparent!(topol[1][1], root.container)
    res   = topol[1][1]

    # 2. Adjust the root position
    state.items[2].t.data = (0.0, -1.0, 0.0)
    ProtoSyn.request_c2i!(pose.state)
    sync!(pose)

    # 1. Define starting positions of the 2 atom template
    x1 = T(0.0)
    y1 = T(0.0)
    x2 = r * cos(π / T(3.0))
    y2 = r * sin(π / T(3.0))

    a0 = Atom("CA011", -1, -1, "C")
    push!(res, a0)
    ProtoSyn.setparent!(a0, root)
    a0s = AtomState()
    a0s.b = 2*y2
    a0s.θ = T(π)
    push!(state, a0s)

    b0 = Atom("CB011", -1, -1, "C")
    push!(res, b0)
    ProtoSyn.setparent!(b0, a0)
    b0s = AtomState()
    b0s.b = r
    b0s.θ = T(-π/3)
    push!(state, b0s)
    

    # # 2. Define movement vectors to copy the 2 atom template
    # # ax, ay = T(2.0) * x2, T(0.0)
    # # bx, by = x2, T(1.5) * r

    # # 3. Populate layer with atoms by copying the original 2 atoms template
    # for j in 1:y
    #     for i in 1:x
    #         println(" i: $i | j: $j")
    #         # dx = i * ax + j * bx
    #         # dy = i * ay + j * by

    #         a0 = Atom("C0$i$j", -1, -1, "C")
    #         ProtoSyn.setparent

    #         atom1_state = AtomState()
    #         atom1_state.b = 
    #         push!(state, atom1_state)
            
    #         a1 = Atom("C$i$(j)1", -1, -1, "C")
    #         ProtoSyn.setparent!(a1, anchor)
    #         push!(a1, res)

    #         atom1_state = AtomState()
    #         atom1_state.t.data = (x2 + dx, y2 + dy, T(0.0))
    #         push!(state, atom1_state)

    #         Atom!(res, "C$i$(j)2", -1, -1, "C")
    #     end
    # end

    topol.id = state.id = genid()
    pose     = Pose(topol, state)
    reindex(pose)
    ProtoSyn.request_i2c!(pose.state)
    sync!(pose)

    return pose
end