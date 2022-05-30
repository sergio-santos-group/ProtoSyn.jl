using ProtoSyn.Units
using LinearAlgebra
using Printf

"""
    generate_carbon_layer(x::Int, y::Int; [r::T = 1.4]) where {T <: AbstractFloat}

The carbon layer is generated using internal coordinates: all atoms are always
placed at a distance `r` and angle of 120° from the parent atom, either with a
clockwise or counter-clockwise rotation in the dihedral angle. For each ring
added, 4 parameters need to be defined:
* The number of atoms to add (N);
* The selected atom to continue X expansion (X);
* The selected atom to continue Y expansion (Y);
* The bonds to add (B).
Such parameters are first compiled in 4 matrixes (1 for each) which are then
consumed to generate the actual atom positions and bonds. This algorithm is,
therefore, sensitive to the order of addition of layers (wither to add Y axis or
X axis first, changes the N, X, Y and B matrixes, as well as rotation of
dihedral angles). For this particular implementation, the X axis is added first,
and all expansion is then added from this layer in the Y axis. As atoms are
added, the next anchor atoms (as selected by the X and Y matrixes) are added to
the A matrix. Each new ring checks this matrix for the next anchor/parent atom.
All dihedral angles on the first atom of each new ring added are rotated by 180°
(default is 0°), except on the first layer where the second atom added is also
rotated by 180°. 6 possible bonding patterns are available:
1. Bond 4th atom to start atom
2. Bond 4th atom to 3rd atom of previous ring
3. Bond 1st atom to 3rd atom of previous ring
4. Bond 4th atom to parent of start atom
5. Perform both 3 and 4 patterns simultaneously
6. Bond 2nd atom to 1st atom of previous ring
Besides this "extra" bonds, all added atoms are bonded to the respective parent
atom. Note that it's guaranteed that the main plane of the generated layer lies
in the x-y plane (z cartesian coordinate = 0.0).
"""
function generate_carbon_layer(x::Int, y::Int; r::T = 1.4) where {T <: AbstractFloat}

    """
    Add an atom to the `parent` atom, named `name`, based on the given internal
    coordinates (distance `b` and dihedral angle `ϕ`, the atomic angle is set to
    180°). If a Residue `residue` isn't explicitly given, use the
    `parent.container`. If `add_bond` is set to `true` (is, by default), add a
    bond record to both `parent` and the newly defined atom.
    """
    function add_atom_from_internal_coordinates!(parent::Atom, name::String, b::T, ϕ::T; residue::Opt{Residue} = nothing, add_bond::Bool = true) where {T <: AbstractFloat}
        a = Atom(name, -1, -1, "C")
        if residue === nothing
            push!(parent.container, a)
        else
            push!(residue, a)
        end
        ProtoSyn.setparent!(a, parent)
        add_bond && ProtoSyn.bond(a, parent)
        as = AtomState()
        as.b = b
        as.θ = 120°
        as.ϕ = ϕ
        push!(state, as)
        return a
    end

    """
    Add "extra" bonds.
    """
    function add_bonds!(added_atoms::Vector{Atom}, previously_added_atoms::Vector{Atom}, start_atom::Atom, pattern::Int)
        if pattern === 0
            nothing
        elseif pattern === 1     # 5 -> start
            ProtoSyn.bond(added_atoms[5], start_atom)
        elseif pattern === 2 # 4 -> 3
            ProtoSyn.bond(added_atoms[4], previously_added_atoms[3])
        elseif pattern === 3 # 1 -> 2
            ProtoSyn.bond(added_atoms[1], previously_added_atoms[2])
        elseif pattern === 4 # 4 -> parent of start
            ProtoSyn.bond(added_atoms[4], start_atom.parent)
        elseif pattern === 5 # pattern 3 + pattern 4
            ProtoSyn.bond(added_atoms[1], previously_added_atoms[2])
            ProtoSyn.bond(added_atoms[4], start_atom.parent)
        elseif pattern === 6 # 2 -> 1
            ProtoSyn.bond(added_atoms[2], previously_added_atoms[1])
        else
            @error "Pattern $pattern not recognized. Use pattern from 0 to 6."
        end
    end

    # 1. Define the pose (all atoms are added to a single residue)
    state    = State()
    topol    = Topology("CRV", -1)
    topol.id = state.id = genid()
    pose     = Pose(topol, state)
    root     = ProtoSyn.root(topol)
    Segment!(topol, "CRV", -1)
    Residue!(topol[1], "CRV", -1)
    ProtoSyn.setparent!(topol[1][1], root.container)

    # 2. Define the correct position for the pose root, so that the plane of the
    # new layer remains in the x-y plane
    x1 = sin(60°)*r
    y1 = cos(60°)*r
    pose.state[root].t.data               = (  -x1,    y1, T(0.0))
    pose.state[root.parent].t.data        = (-2*x1, T(0.0), T(0.0))
    pose.state[root.parent.parent].t.data = (-3*x1,    y1, T(0.0))
    ProtoSyn.request_c2i!(pose.state, all = true)
    sync!(pose)

    # 3. Define the starting atom (conneted to root)
    a0 = add_atom_from_internal_coordinates!(root, "C0", r, 180°,
        residue = topol[1, 1], add_bond = false)

    # 4. Generate N matrix (Number of addition) & generate the B matrix (bonding
    # patterns on each ring)
    N = zeros(Int, y, x) .+ 2
    N[y, :] .+= 2
    N[y, 1]  += 1

    B = zeros(Int, y, x)
    B[y, :] .+= 1
    B[y, 2:end] .+= 1

    f = y % 2 === 0 ? (i) -> i % 2 === 0 : (i) -> i % 2 !== 0
    for i in 1:(y-1)
        p1 = f(i) ? [i, 1] : [i, x]
        N[p1...] += 2
        B[p1...] += f(i) ? 4 : 5
        p2 = f(i) ? [i, 2:x] : [i, (2:(x-1))]
        B[p2...] .+= f(i) ? 6 : 3
    end

    # 5. Generate X matrix (X axis expansion atom selection)
    X = zeros(Int, y, x)
    X[y, 1:(x-1)] .+= 2

    # 6. Generate Y matrix (Y axis expansion atom selection)
    Y = zeros(Int, y, x)
    Y[2:y, :] .+= 2
    Y[y, :] .+= 2

    # 7. Generate the A matrix (starting atoms), initialized by the a0
    a = mapreduce(permutedims, vcat, [[nothing for _ in 1:x] for _ in 1:y])
    A = Matrix{Opt{Atom}}(a)
    A[y, 1] = a0

    # 8. Consume the defined matrixes
    previously_added_atoms = Vector{Atom}([])
    for i in y:-1:1
        for j in 1:x
            @debug "($i, $j) => N: $(N[i, j]) | X: $(X[i, j]) | Y: $(Y[i, j])"
            start_atom  = A[i, j]
            a           = A[i, j]
            added_atoms = Vector{Atom}([])
            @debug " Start atom: $start_atom"
            for n in 1:N[i, j]
                if i === y
                    ϕ = (n in [1, 2]) ? 180° : 0°
                else
                    ϕ = n !== 1 ? 0° : 180°
                end
                @debug @printf(" Adding atom %d (Name: %s | Rotation %6.2f)\n", n, "C", rad2deg(ϕ))
                a = add_atom_from_internal_coordinates!(a, "C", r, ϕ)
                push!(added_atoms, a)
            end
            
            @debug "  Adding next X-axis atom: $(X[i, j] !== 0)"
            if X[i, j] !== 0
                @debug "  Next X atom: $(added_atoms[X[i, j]])"
                A[i, j+1] = added_atoms[X[i, j]]
            end

            @debug "  Adding next Y-axis atom: $(Y[i, j] !== 0 && i-1 > 0)"
            if Y[i, j] !== 0 && i-1 > 0
                @debug "  Next Y atom: $(added_atoms[Y[i, j]])"
                A[i-1, j] = added_atoms[Y[i, j]]
            end

            pattern = B[i, j]
            @debug "Pattern: $pattern"
            add_bonds!(added_atoms, previously_added_atoms, start_atom, pattern)
            previously_added_atoms = copy(added_atoms)
        end
    end

    # 9. Reindex pose and sync! internal to cartesian coordinates
    reindex(pose)
    ProtoSyn.request_i2c!(pose.state)
    sync!(pose)

    return pose
end


"""
# TODO: Documentation
"""
function generate_carbon(x::Int, y::Int, z::Int = 0; r::T = 1.4, d::T = 3.4) where {T <: AbstractFloat}
    pose = generate_carbon_layer(x, y, r = r)
    _v1 = pose.state[pose.graph[1, 1, 1]].t - pose.state[pose.graph[1, 1, 2]].t
    _v2 = pose.state[pose.graph[1, 1, 6]].t - pose.state[pose.graph[1, 1, 1]].t
    v1 = normalize(cross(_v1, _v2))
    for i in 2:z
        layer = generate_carbon_layer(x, y, r = r)
        
        N = layer.state.size
        D = (i-1)*d
        println("D: $D")
        println("Before: $(layer.state[1].t[3])")
        layer.state.x[:, 1:N] = layer.state.x[:, 1:N] .+ (v1 .* [D, D, D])
        
        ProtoSyn.request_c2i!(layer.state)
        sync!(layer)
        println("After: $(layer.state[1].t[3])")



        if i % 2 === 0 # set multiple layer offset
            println("OFFSET")
            v2 = layer.state[layer.graph[1, 1, 1]].t .- ProtoSyn.center_of_mass(layer, aid"1:6")[:]
            layer.state.x[:, 1:N] = layer.state.x[:, 1:N] .+ (v2 .* [1, 1, 0])
        end

        ProtoSyn.request_c2i!(layer.state)
        sync!(layer)

        ProtoSyn.merge!(pose, layer)
    end

    return pose
end


"""
# TODO: Documentation
"""
function generate_porosity(pose::Pose, pore_fraction::T, clean_sweeps::Int = 15; random::Bool = false) where {T <: AbstractFloat}

    for atom in eachatom(pose.graph)
        as = pose.state[atom]
        if random
            random_scale = 0.1
            rs = ((rand()*2)-1) * random_scale
            noise = perlinoctaves(T(as.t[1]), T(as.t[2]), T(as.t[3]), rand([2, 3, 4, 5]), ((rand()*2)+1) * random_scale) + rs
        else
            noise = perlin(T(as.t[1]), T(as.t[2]), T(as.t[3]))
        end
        
        if noise > pore_fraction
            ProtoSyn.pop_atom!(pose, atom, keep_downstream_position = false)
        end
    end

    for _ in 1:clean_sweeps
        for atom in eachatom(pose.graph)
            if length(atom.bonds) <= 1
                ProtoSyn.pop_atom!(pose, atom, keep_downstream_position = false)
            end
        end
    end

    pose.state.i2c = false
    ProtoSyn.request_c2i!(pose.state)
    sync!(pose)
    return pose
end