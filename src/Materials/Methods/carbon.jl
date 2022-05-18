function new_layer(x::Int, y::Int, r::T = 1.4) where {T <: AbstractFloat}

    state = State()
    topol = Topology("CRV", -1)
    Segment!(topol, "CRV", -1)
    Residue!(topol[1], "CRV", -1)
    res   = topol[1][1]

    # 1. Define starting positions of the 2 atom template
    x1 = T(0.0)
    y1 = r
    x2 = r * cos(π / T(6.0))
    y2 = r * sin(π / T(6.0))

    # 2. Define movement vectors to copy the 2 atom template
    ax, ay = T(2.0) * x2, T(0.0)
    bx, by = x2, T(1.5) * r

    # 3. Populate layer with atoms by copying the original 2 atoms template
    for j in 1:(y+1)
        for (n, i) in enumerate((-j ÷ 2):(x - j ÷ 2))
            dx = i * ax + j * bx
            dy = i * ay + j * by
            if !(y % 2 == 0) && (j == y) && (n == 0)
                nothing
            else
                atom1_state = AtomState()
                atom1_state.t.data = (x1 + dx, y1 + dy, T(0.0))
                push!(state, atom1_state)
                
                Atom!(res, "C", -1, -1, "C")
            end
            atom1_state = AtomState()
            atom1_state.t.data = (x2 + dx, y2 + dy, T(0.0))
            push!(state, atom1_state)

            Atom!(res, "C", -1, -1, "C")
        end
    end

    topol.id = state.id = genid()
    pose     = Pose(topol, state)
    root     = ProtoSyn.root(pose.graph)
    for atom in eachatom(pose.graph)
        ProtoSyn.setparent!(atom, root)
    end
    reindex(pose)
    ProtoSyn.request_c2i!(pose.state, all = true)
    sync!(pose)

    return pose
end