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
therefore, sensitive to the order of addition of layers (whether to add Y axis
or X axis first, changes the N, X, Y and B matrixes, as well as rotation of
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

# See also
[`generate_carbon`](@ref)

# Examples
```
julia> pose = ProtoSyn.Materials.generate_carbon_layer(5, 5)
Pose{Topology}(Topology{/CRV:65069}, State{Float64}:
 Size: 70
 i2c: false | c2i: false
 Energy: Dict(:Total => Inf)
)
```
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
    y1 = sin(60°)*r
    x1 = cos(60°)*r
    pose.state[root].t.data               = (    -x1,    -y1, T(0.0))
    pose.state[root.parent].t.data        = (  -x1-r,    -y1, T(0.0))
    pose.state[root.parent.parent].t.data = (-2*x1-r, T(0.0), T(0.0))
    root.ascendents = (0, -1, -2, -3)
    root.parent.ascendents = (-1, -2, -3, -4)
    root.parent.parent.ascendents = (-2, -3, -4, -5)
    ProtoSyn.request_c2i!(pose.state, all = true)
    sync!(pose)
    
    # 3. Define the starting atom (conneted to root)
    a0 = add_atom_from_internal_coordinates!(root, "C0", r, 0°,
        residue = topol[1, 1], add_bond = false)

    reindex(pose.graph, set_ascendents = true)
    reindex(pose.state)
    ProtoSyn.request_i2c!(pose.state, all = true)
    sync!(pose)
    ProtoSyn.request_c2i!(pose.state, all = true)
    sync!(pose)

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
            @info "($i, $j) => N: $(N[i, j]) | X: $(X[i, j]) | Y: $(Y[i, j])"
            start_atom  = A[i, j]
            a           = A[i, j]
            added_atoms = Vector{Atom}([])
            @info " Start atom: $start_atom"
            for n in 1:N[i, j]
                if i === y
                    ϕ = (n in [1, 2]) ? 180° : 0°
                else
                    ϕ = n !== 1 ? 0° : 180°
                end
                @info @printf(" Adding atom %d (Name: %s | Rotation %6.2f)\n", n, "C", rad2deg(ϕ))
                a = add_atom_from_internal_coordinates!(a, "C", r, ϕ)
                push!(added_atoms, a)
            end
            
            @info "  Adding next X-axis atom: $(X[i, j] !== 0)"
            if X[i, j] !== 0
                @info "  Next X atom: $(added_atoms[X[i, j]])"
                A[i, j+1] = added_atoms[X[i, j]]
            end

            @info "  Adding next Y-axis atom: $(Y[i, j] !== 0 && i-1 > 0)"
            if Y[i, j] !== 0 && i-1 > 0
                @info "  Next Y atom: $(added_atoms[Y[i, j]])"
                A[i-1, j] = added_atoms[Y[i, j]]
            end

            # Add non-parenthood bonds
            pattern = B[i, j]
            @info "Pattern: $pattern"
            add_bonds!(added_atoms, previously_added_atoms, start_atom, pattern)
            previously_added_atoms = copy(added_atoms)
        end
    end

    # 9. Reindex pose and sync! internal to cartesian coordinates
    reindex(pose)
    ProtoSyn.request_i2c!(pose.state, all = true)
    sync!(pose)

    return pose
end


"""
    generate_carbon(x::Int, y::Int, [z::Int = 0]; [r::T = 1.4], [d::T = 3.4]) where {T <: AbstractFloat}

Generate a carbon microcrystallite with multiple carbon sheet layers (built
using [`generate_carbon_layer`](@ref) with `x` and `y` dimensions for the number
of carbon rings in the x- and y-axis, respectivelly, and `r` distance between
atoms, in Å). The generated layers are copied `z` times (in the x-axis) with a
distance of `d` Å. Every odd layer is displaced to match π-π stacking chemistry
in carbon sheets.

# See also
[`generate_porosity`](@ref) [`generate_carbon_from_file`](@ref)

# Examples
```
julia> pose = ProtoSyn.Materials.generate_carbon(5, 5, 5)
Pose{Topology}(Topology{/CRV:23442}, State{Float64}:
 Size: 350
 i2c: false | c2i: false
 Energy: Dict(:Total => Inf)
)
```
"""
function generate_carbon(x::Int, y::Int, z::Int = 0; r::T = 1.4, d::T = 3.4) where {T <: AbstractFloat}
    pose = copy(generate_carbon_layer(x, y, r = r))
    ProtoSyn.request_c2i!(pose.state, all = false) # ! Not sure why this is necessary.
    sync!(pose)

    v1 = [0.0, 0.0, 1.0]

    for i in 2:z
        layer = copy(generate_carbon_layer(x, y, r = r))
        
        N = layer.state.size
        D = (i-1)*d
        layer.state.x[:, 1:N] = layer.state.x[:, 1:N] .+ (v1 .* [D, D, D])
        

        if i % 2 === 0 # set multiple layer offset
            v2 = [r/2, cos(30°)*r, 0.0]
            layer.state.x[:, 1:N] = layer.state.x[:, 1:N] .+ (v2 .* [1, 1, 0])
        end

        ProtoSyn.request_c2i!(layer.state, all = false)
        sync!(layer)

        ProtoSyn.merge!(pose, layer)
    end

    ProtoSyn.request_i2c!(pose.state, all = true)
    sync!(pose)

    return pose
end


"""
    generate_porosity(pose::Pose, pore_fraction::T; [clean_sweeps::Int = 15], [random::Bool = false], [neat_indexation::Bool = false]) where {T <: AbstractFloat}

Generate porosity on the given [`Pose`](@ref) `pose` by removing [`Atom`](@ref)
instances according to a generated Perlin noise (see the [`perlin`](@ref)
method). The `pore_fraction` field ajusts the amount of [`Atom`](@ref) instances
removed, as a value between 0.0 and 1.0.

!!! ukw "Note:"
    The `pose_fraction` field cannot be roughly translated as a direct percentage value of [`Atom`](@ref) instances to remove. Instead, this value refers to the [`perlin`](@ref) noise level to consider for atom removal, and as such scales exponentially. As a baseline, a `pore_fraction` of 0.4 removes approximately 50% of the [`Atom`](@ref) instances in the given [`Pose`](@ref). 

For high `pore_fraction` values, the [`Atom`](@ref) removal may leave
unconnected chains of [`Atom`](@ref) instances in each [`Segment`](@ref). Such
cases should be manually verified. Since this function is intended to be applied
to carbon microcrystallites, these hanging atoms may constitute carbon atoms
with wrong valency numbers. This method automatically performs `clean_sweeps`
and removes any atom with 1 or 0 bonds. If `random` is set to `true` (`false`,
by default), the generated noise is randomized. If `neat_indexation` is set to
`true` (`false`, by default), the
[`sort_atoms_by_graph!`](@ref ProtoSyn.sort_atoms_by_graph!) method is applied
to reindex [`Atom`](@ref) instances and sort them according to the new graph
(after removing [`Atom`](@ref) instances).

# See also
[`generate_carbon`](@ref)

# Examples
```
julia> ProtoSyn.Materials.generate_porosity(pose, 0.405)
Pose{Topology}(Topology{/CRV:36031}, State{Float64}:
 Size: 190
 i2c: false | c2i: false
 Energy: Dict(:Total => Inf)
)
```
"""
function generate_porosity(pose::Pose, pore_fraction::T; clean_sweeps::Int = 15, random::Bool = false, neat_indexation::Bool = false) where {T <: AbstractFloat}

    for atom in eachatom(pose.graph)
        as = pose.state[atom]
        if random
            random_scale = 0.1
            rs = ((rand()*2)-1) * random_scale
            noise = perlinoctaves(T(as.t[1]), T(as.t[2]), T(as.t[3]), rand([2, 3, 4, 5]), ((rand()*2)+1) * random_scale) + rs
        else
            noise = perlin(T(as.t[1]), T(as.t[2]), T(as.t[3]))
        end
        
        if noise > (T(1.0) - pore_fraction)
            @info " Pore generation in carbon sheet: Popping atom $atom"
            ProtoSyn.pop_atom!(pose, atom, keep_downstream_position = true)
        end
    end

    # Clean sweeps deletes "hanging" atoms (with only 1 or no bonds). This is
    # performed N times, as from the previous attempt new hanging atoms might
    # have been generated.
    for _ in 1:clean_sweeps
        for atom in eachatom(pose.graph)
            if length(atom.bonds) <= 1
                ProtoSyn.pop_atom!(pose, atom, keep_downstream_position = true)
            end
        end
    end

    if ProtoSyn.count_atoms(pose.graph) === 0
        @error "Porosity too high! ProtoSyn removed all atoms during porosity generation!"
    end

    # During pop_atom!, multiple atoms were set to be child of root. A new
    # infer_parenthood! round makes it so that all atoms are children of other
    # atoms (only N root children remain, for the N segments in the pose, as
    # long as all segments are contiguous. Non-contiguous sub-segments will
    # still be children of root and may cause problems downstream).
    for segment in eachsegment(pose.graph)
        ProtoSyn.count_atoms(segment) === 0 && continue
        ProtoSyn.infer_parenthood!(pose.graph, overwrite = true, start = segment[1, 1], linear_aromatics = false) 
    end

    # Since the graph structure has changes, a new update of the internal
    # coordinates needs to be performed
    ProtoSyn.request_c2i!(pose.state, all = true)
    sync!(pose)

    # Given the parenhood remixes performed, the current indexation might have
    # been lost. If `neat_indexation` is set to true, the `sort_atoms_by_graph!`
    # attempts to recover a more "linear" indexation.
    if neat_indexation
        for segment in eachsegment(pose.graph)
            ProtoSyn.count_atoms(segment) === 0 && continue
            ign_sele = !SerialSelection{Segment}(segment.id, :id)
            ProtoSyn.sort_atoms_by_graph!(pose, start = segment[1, 1], search_algorithm = ProtoSyn.BFS, ignore_selection = ign_sele)
        end

        reindex(pose.graph, set_ascendents = true)
        reindex(pose.state)

        ProtoSyn.request_c2i!(pose.state, all = true)
        sync!(pose)
    else
        reindex(pose.graph, set_ascendents = true)
        reindex(pose.state)
    end

    return pose
end


"""
    generate_carbon_from_file(filename::String, output::Opt{String} = nothing)

Generate a functionalized carbon model from a .yml input file (named
`filename`). If `output` is set to a `String` (`nothing`, by default), outputs
the generated carbon model to a file. For an input file example, check the
`resources/Materials/carbon.yml` file.

# See also
[`functionalize!`](@ref)

# Examples
```
julia> pose = ProtoSyn.Materials.generate_carbon_from_file("carbon.yml")
Pose{Topology}(Topology{/CRV:40141}, State{Float64}:
 Size: 960
 i2c: false | c2i: false
 Energy: Dict(:Total => Inf)
)

```
"""
function generate_carbon_from_file(filename::String, output::Opt{String} = nothing)

    # 1. Read file
    data = ProtoSyn.read_yml(filename)

    # 2. Generate carbon
    shape = data["Shape"]
    pose  = generate_carbon(shape["x"], shape["y"], shape["z"])

    # 3. Generate porosity
    if "porosity" in keys(shape) && shape["porosity"] > 0.0
        generate_porosity(pose, shape["porosity"], random = true, neat_indexation = true)
    end

    # 4. Add functional groups
    if "Functional-groups" in keys(data)
        f = Dict{Fragment, eltype(pose.state)}()
        for (fcn, value) in data["Functional-groups"]
            f[ProtoSyn.getvar(ProtoSyn.modification_grammar, fcn)] = value
        end
        functionalize!(pose, f)
    end

    # 5. Add hydrogens
    if "Hydrogens" in keys(data) && "saturate" in keys(data["Hydrogens"])
        if data["Hydrogens"]["saturate"]
            add_hydrogens!(pose, ProtoSyn.modification_grammar, nothing)
        end
    end

    # 5. Write output to file
    if output !== nothing
        println("Outputing results to $output")
        ProtoSyn.write(pose, output)
    end

    println("All tasks done!")

    return pose
end