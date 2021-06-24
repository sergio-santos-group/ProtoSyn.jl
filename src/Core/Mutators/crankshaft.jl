using Printf

"""
    CrankshaftMutator(angle_sampler::Function, p_mut::AbstractFloat, step_size::AbstractFloat, selection::Opt{AbstractSelection}, inc_last_res::Opt{AbstractSelection})

Return a [`CrankshaftMutator`](@ref) instance. This `AbstractMutator` is a
_functor_, called with the following signature:

```
(crankshaft_mutator::CrankshaftMutator)(pose::Pose)
```

The [`CrankshaftMutator`](@ref) `AbstractMutator` loops through all pairs of two
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
`AbstractSelection` `selection` is provided, only the selected [`Atom`](@ref)
instances are looped over. If the given `AbstractSelection` `selection` is not
of selection type [`Atom`](@ref), it will be promoted to this type (using
[`promote`](@ref ProtoSyn.promote) with default aggregator `any`). If an
`AbstractSelection` `inc_last_res` is provided, all atoms in this
`AbstractSelection` and in the last residue considered for the crankshaft
movement will be included in the rotation. As an example, for a correct
crankshaft movement in peptidic structures, the sidechain of the last residue
should also be included in the rotation. In this example, `inc_last_res` should
be `!(an"^CA\$|^N\$|^C\$|^H\$|^O\$"r)`. Note that the
[`CrankshaftMutator`](@ref) syncs any pending internal to cartesian coordinate
conversion (using the [`i2c!`](@ref ProtoSyn.i2c!) method). Requests cartesian
to internal coordinates conversion (using
[`request_c2i!`](@ref ProtoSyn.request_c2i!) method). Does not [`sync!`](@ref)
the given [`Pose`](@ref) afterwards.

The [`CrankshaftMutator`](@ref) `AbstractMutator` can also be optionally called
using the following signature, in which case only the provided list of
[`Atom`](@ref) instances will be considered for the application of this
`AbstractMutator`.

```
(crankshaft_mutator::CrankshaftMutator)(pose::Pose, atoms::Vector{Atom})
```

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
⚯  Crankshaft Mutator:
+----------------------------------------------------------------------+
| Index | Field                       | Value                          |
+----------------------------------------------------------------------+
| 1     | angle_sampler               | Function randn                 |
| 2     | p_mut                       | 1.0000                         |
| 3     | step_size                   | 1.0000                         |
+----------------------------------------------------------------------+
 ○  Selection: Not Set
 ○  Include on last residue: Not Set

julia> ProtoSyn.Mutators.CrankshaftMutator(randn, 0.01, 1.0, an"CA", !(an"^CA\$|^N\$|^C\$|^H\$|^O\$"r))
⚯  Crankshaft Mutator:
+----------------------------------------------------------------------+
| Index | Field                       | Value                          |
+----------------------------------------------------------------------+
| 1     | angle_sampler               | Function randn                 |
| 2     | p_mut                       | 0.0100                         |
| 3     | step_size                   | 1.0000                         |
+----------------------------------------------------------------------+
 ● Selection: Set
 └── FieldSelection › Atom.name = CA
 ● Include on last residue: Set
 └── UnarySelection ❯ ! "not" (Atom)
      └── FieldSelection › Atom.name = r"^CA\$|^N\$|^C\$|^H\$|^O\$"
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
        if ProtoSyn.selection_type(crankshaft_mutator.selection) !== Atom
            sele  = ProtoSyn.promote(crankshaft_mutator.selection, Atom) # default aggregator is `any`
            atoms = sele(pose, gather = true)
        else
            atoms = crankshaft_mutator.selection(pose, gather = true)
        end
    end
    
    crankshaft_mutator(pose, atoms)
end


function (crankshaft_mutator::CrankshaftMutator)(pose::Pose, atoms::Vector{Atom})

    # CrankshaftMutator requires updated cartesian coordinates
    ProtoSyn.i2c!(pose.state, pose.graph) # Checks pose.state.i2c flag inside

    for (i, atom_i) in enumerate(atoms[1:(end-1)])
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

function Base.show(io::IO, cm::CrankshaftMutator, level_code::Opt{LevelCode} = nothing)
    level_code = level_code === nothing ? LevelCode() : level_code
    lead       = ProtoSyn.get_lead(level_code)
    inner_lead = ProtoSyn.get_inner_lead(level_code)

    println(io, lead*"⚯  Crankshaft Mutator:")
    println(io, inner_lead*"+"*repeat("-", 70)*"+")
    @printf(io, "%s| %-5s | %-27s | %-30s |\n", inner_lead, "Index", "Field", "Value")
    println(io, inner_lead*"+"*repeat("-", 70)*"+")
    @printf(io, "%s| %-5d | %-27s | %-30s |\n", inner_lead, 1, "angle_sampler", "Function $(cm.angle_sampler)")
    @printf(io, "%s| %-5d | %-27s | %-30.4f |\n", inner_lead, 2, "p_mut", cm.p_mut)
    @printf(io, "%s| %-5d | %-27s | %-30.4f |\n", inner_lead, 3, "step_size", cm.step_size)
    println(io, inner_lead*"+"*repeat("-", 70)*"+")
    
    if cm.selection !== nothing
        println(io, inner_lead*" ● Selection: Set")
        Base.show(io, cm.selection, vcat(level_code, 4))
    else
        println(io, inner_lead*" ○  Selection: Not Set")
    end

    if cm.inc_last_res !== nothing
        println(io, inner_lead*" ● Include on last residue: Set")
        Base.show(io, cm.inc_last_res, vcat(level_code, 4))
    else
        println(io, inner_lead*" ○  Include on last residue: Not Set")
    end
end