using Printf

"""
    CrankshaftMutator(angle_sampler::Function, p_mut::AbstractFloat, step_size::AbstractFloat, selection::Opt{AbstractSelection}, inc_last_res::Opt{AbstractSelection})

Return a [`CrankshaftMutator`](@ref) instance. This `AbstractMutator` is a
_functor_, called with the following signature:

```
(crankshaft_mutator::CrankshaftMutator)(pose::Pose)
```

The [`CrankshaftMutator`](@ref) `AbstractMutator` loop through all pairs of two
atoms possible, in the given [`Pose`](@ref). For each pair, if a random number
(`rand()`) is bellow a given probability of mutation `p_mut`, both atoms are
selected for a crankshaft movement (therefore a higher `p_mut` value applies a
larger number of crankshaft movements per call). A crankshaft movement is a
concerted rotation around a virtual axis formed between the two selected atoms,
by a certain angle given by the `angle_sampler` function and multiplied by the
`step_size` value. The rotation is applied to the region between the two atoms,
by travelling the [`Pose`](@ref) [Graph](@ref graph-types) (See
[`travel_graph`](@ref)), by calculating and applying a rotation matrix (See
[`rotation_matrix_from_axis_angle`](@ref ProtoSyn.rotation_matrix_from_axis_angle)).
Note that a new angle for the crankshaft movement is sampled for each selected
pair of [`Atom`](@ref) instances. If an `AbstractSelection` `selection` is
provided, only the selected [`Atom`](@ref) instances are looped over. If an
`AbstractSelection` `inc_last_res` is provided, all atoms in this
`AbstractSelection` and in the last residue considered for the crankshaft
movement will be included in the rotation. As an example, for a correct
crankshaft movement in peptidic structures, the sidechain of the last residue
should also be included in the rotation. In this example, `inc_last_res` should
be `!(an"^CA\$|^N\$|^C\$|^H\$|^O\$"r)`. Requests cartesian to internal
coordinates conversion (using [`request_c2i!`](@ref ProtoSyn.request_c2i!)
method). Does not [`sync!`](@ref) the given [`Pose`](@ref).

# Fields
* `angle_sampler::Function` - Should return a `Float` angle value (in radians). Is called with no input arguments;
* `p_mut::AbtractFloat` - Compared against a `rand()` call, applies this Mutator to [`Atom`](@ref) instance pairs where `rand() < p_mut`;
* `step_size::AbstractFloat` - Multiplies the sampled angle by this value;
* `selection::Opt{AbstractSelection}` - If given, this Mutator will only loop over the selected [`Atom`](@ref) instances;
* `inc_last_res::Opt{AbstractSelection}` - If given, all atoms in this Mutator and in the last residue considered for the crankshaft movement will be included in the rotation.

# See also
[`DihedralMutator`](@ref)

# Examples
```jldoctest
julia> ProtoSyn.Mutators.CrankshaftMutator(randn, 1.0, 1.0, nothing, nothing)
  Crankshaft:
+--------------------------------------------------------------------+
| Index | Field                     | Value                          |
+--------------------------------------------------------------------+
| 1     | angle_sampler             | Function randn                 |
| 2     | p_mut                     | 1.0000                         |
| 3     | step_size                 | 1.0000                         |
| 4     | selection                 | Not set                        |
| 4     | inc_last_res              | Not set                        |
+--------------------------------------------------------------------+

julia> ProtoSyn.Mutators.CrankshaftMutator(randn, 0.01, 1.0, an"CA", !(an"^CA\$|^N\$|^C\$|^H\$|^O\$"r))
  Crankshaft:
+--------------------------------------------------------------------+
| Index | Field                     | Value                          |
+--------------------------------------------------------------------+
| 1     | angle_sampler             | Function randn                 |
| 2     | p_mut                     | 0.0100                         |
| 3     | step_size                 | 1.0000                         |
| 4     | selection                 | Set: FieldSelection            |
| 4     | inc_last_res              | Set: FieldSelection            |
+--------------------------------------------------------------------+
```
"""
mutable struct CrankshaftMutator <: AbstractMutator
    angle_sampler::Function # Should return a float in radians
    p_mut::AbstractFloat
    step_size::AbstractFloat
    selection::Opt{AbstractSelection}
    inc_last_res::Opt{AbstractSelection} # !(an"^CA$|^N$|^C$|^H$|^O$"r) in peptides
    # Note: the "inc_last_res" parameter is a selection. All atoms in
    # the selection AND in the last residue considered for the crankshaft
    # movement WILL BE included in the rotation.
end

function (crankshaft_mutator::CrankshaftMutator)(pose::Pose)
    if crankshaft_mutator.selection === nothing
        atoms = collect(eachatom(pose.graph))
    else
        atoms = crankshaft_mutator.selection(pose, gather = true)
    end
    
    crankshaft_mutator(pose, atoms)
end


function (crankshaft_mutator::CrankshaftMutator)(pose::Pose, atoms::Vector{Atom})

    for (i, atom_i) in enumerate(atoms)
        for atom_j in atoms[(i+1):end]
            if rand() < crankshaft_mutator.p_mut
                # 1) Get angle
                ∠      = crankshaft_mutator.angle_sampler()
                ∠     *= crankshaft_mutator.step_size

                # 2) Get axis
                # Note 1: Ai is the pivot
                pivot  = pose.state[atom_i].t
                axis   = collect(pose.state[atom_j].t - pivot)

                # 3) Get rotation matrix from axis and angle
                rmat   = ProtoSyn.rotation_matrix_from_axis_angle(axis, ∠)

                # 4) Apply the rotation matrix to all atoms between Ai and Aj
                # Note 2: inc_last_res selection optionally adds a subset of
                # atom from the last residue in the rotation. In the case of
                # proteins, for example, it should be the sidechain, as in
                # the selection !(an"^CA$|^N$|^C$|^H$|^O$"r)
                ids    = ProtoSyn.ids(travel_graph(atom_i, atom_j))
                atoms  = crankshaft_mutator.inc_last_res(atom_j.container, gather = true)
                ids2   = ProtoSyn.ids(atoms)
                idxs   = vcat(ids, ids2)
                M      = pose.state.x[:, idxs]
                result = (rmat * (M .- pivot)) .+ pivot
                pose.state.x[:, idxs] = result
                ProtoSyn.request_c2i!(pose.state, all = true)
            end
        end
    end
end

function Base.show(io::IO, cm::CrankshaftMutator)
    println("  Crankshaft:")
    println(io, "+"*repeat("-", 68)*"+")
    @printf(io, "| %-5s | %-25s | %-30s |\n", "Index", "Field", "Value")
    println(io, "+"*repeat("-", 68)*"+")
    @printf(io, "| %-5d | %-25s | %-30s |\n", 1, "angle_sampler", "Function $(cm.angle_sampler)")
    @printf(io, "| %-5d | %-25s | %-30.4f |\n", 2, "p_mut", cm.p_mut)
    @printf(io, "| %-5d | %-25s | %-30.4f |\n", 3, "step_size", cm.step_size)
    if cm.selection === nothing
        @printf(io, "| %-5d | %-25s | %-30s |\n", 4, "selection", "Not set")
    else
        @printf(io, "| %-5d | %-25s | %-30s |\n", 4, "selection", "Set: $(typeof(cm.selection).name.name)")
    end
    if cm.inc_last_res === nothing
        @printf(io, "| %-5d | %-25s | %-30s |\n", 4, "inc_last_res", "Not set")
    else
        @printf(io, "| %-5d | %-25s | %-30s |\n", 4, "inc_last_res", "Set: $(typeof(cm.selection).name.name)")
    end
    println(io, "+"*repeat("-", 68)*"+")
end