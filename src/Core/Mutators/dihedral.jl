using Printf

"""
    DihedralMutator(angle_sampler::Function, p_mut::AbstractFloat, step_size::AbstractFloat, selection::Opt{AbstractSelection})

Return a new [`DihedralMutator`](@ref) instance. This `AbstractMutator` is a
_functor_, called with the following signature:

```
(dihedral_mutator::DihedralMutator)(pose::Pose)
```

The [`DihedralMutator`](@ref) `AbstractMutator` loops through all the
[`Atom`](@ref) instances in a given [`Pose`](@ref) and applies a new dihedral
angle if a random number (`rand()`) is bellow a given probability of mutation
`p_mut` (therefore a higher `p_mut` value applies a larger number of dihedral
rotations per call). The applied rotation is an angle value (in radians),
sampled by calling `angle_sampler` and multiplied by the `step_size` value. The
resulting value is then added to the selected [`Atom`](@ref)`.Δϕ`. Note that a
new dihedral rotation is sampled for each selected [`Atom`](@ref) instance. If
an `AbstractSelection` `selection` is provided, only the selected [`Atom`](@ref)
instances are looped over. If the given `AbstractSelection` `selection` is not
of selection type [`Atom`](@ref), it will be promoted to this type (using
[`promote`](@ref ProtoSyn.promote) with default aggregator `any`). Note that the
[`DihedralMutator`](@ref) syncs any pending cartesian to internal coordinate
conversion (using the [`c2i!`](@ref ProtoSyn.c2i!) method). Requests internal to
cartesian coordinates conversion (using
[`request_i2c!`](@ref ProtoSyn.request_i2c!) method). Does not [`sync!`](@ref)
the given [`Pose`](@ref) afterwards.

The [`DihedralMutator`](@ref) `AbstractMutator` can also be optionally called
using the following signature, in which case only the provided list of
[`Atom`](@ref) instances will be considered for the application of this
`AbstractMutator`.

```
(dihedral_mutator::DihedralMutator)(pose::Pose, atoms::Vector{Atom})
```

# Fields

* `angle_sampler::Function` - Should return a `Float` angle value (in radians). Is called with no input arguments;
* `p_mut::AbtractFloat` - Compared against a `rand()` call, applies this Mutator to [`Atom`](@ref) instances where `rand() < p_mut`;
* `step_size::AbstractFloat` - Multiplies the sampled angle by this value;
* `selection::Opt{AbstractSelection}` - If given, this Mutator will only loop over the selected [`Atom`](@ref) instances.

# See also
[`CrankshaftMutator`](@ref)

# Examples
```jldoctest
julia> ProtoSyn.Mutators.DihedralMutator(randn, 1.0, 1.0, nothing)
⚯  Dihedral Mutator:
+----------------------------------------------------------------------+
| Index | Field                       | Value                          |
+----------------------------------------------------------------------+
| 1     | angle_sampler               | Function randn                 |
| 2     | p_mut                       | 1.0000                         |
| 3     | step_size                   | 1.0000                         |
+----------------------------------------------------------------------+
 ○  Selection: Not Set

julia> ProtoSyn.Mutators.DihedralMutator(randn, 0.05, 1.0, an"CA\$|C\$"r)
⚯  Dihedral Mutator:
+----------------------------------------------------------------------+
| Index | Field                       | Value                          |
+----------------------------------------------------------------------+
| 1     | angle_sampler               | Function randn                 |
| 2     | p_mut                       | 0.0500                         |
| 3     | step_size                   | 1.0000                         |
+----------------------------------------------------------------------+
 ● Selection: Set
 └── FieldSelection › Atom.name = r"CA\$|C\$"
```
"""
mutable struct DihedralMutator <: AbstractMutator
    angle_sampler::Function # Should return a float in radians
    p_mut::AbstractFloat
    step_size::AbstractFloat
    selection::Opt{AbstractSelection}
end

function (dihedral_mutator::DihedralMutator)(pose::Pose)
    if dihedral_mutator.selection === nothing
        atoms = collect(eachatom(pose.graph))
    else
        if ProtoSyn.selection_type(dihedral_mutator.selection) !== Atom
            sele  = ProtoSyn.promote(dihedral_mutator.selection, Atom) # default aggregator is `any`
            atoms = sele(pose, gather = true)
        else
            atoms = dihedral_mutator.selection(pose, gather = true)
        end
    end
    
    dihedral_mutator(pose, atoms)
end

function (dihedral_mutator::DihedralMutator)(pose::Pose, atoms::Vector{Atom})
    
    # DihedralMutator requires updated internal coordinates
    sync!(pose)
    # ProtoSyn.c2i!(pose.state, pose.graph) # Checks pose.state.c2i flag inside

    for atom in atoms
        if rand() < dihedral_mutator.p_mut
            ∠ = dihedral_mutator.angle_sampler() * dihedral_mutator.step_size
            ProtoSyn.rotate_dihedral!(pose.state, atom, ∠) # Already requests i2c
        end
    end
end

function Base.show(io::IO, dm::DihedralMutator, level_code::Opt{LevelCode} = nothing)
    level_code = level_code === nothing ? LevelCode() : level_code
    lead       = ProtoSyn.get_lead(level_code)
    inner_lead = ProtoSyn.get_inner_lead(level_code)

    println(io, lead*"⚯  Dihedral Mutator:")
    println(io, inner_lead*"+"*repeat("-", 70)*"+")
    @printf(io, "%s| %-5s | %-27s | %-30s |\n", inner_lead, "Index", "Field", "Value")
    println(io, inner_lead*"+"*repeat("-", 70)*"+")
    @printf(io, "%s| %-5d | %-27s | %-30s |\n", inner_lead, 1, "angle_sampler", "Function $(dm.angle_sampler)")
    @printf(io, "%s| %-5d | %-27s | %-30.4f |\n", inner_lead, 2, "p_mut", dm.p_mut)
    @printf(io, "%s| %-5d | %-27s | %-30.4f |\n", inner_lead, 3, "step_size", dm.step_size)
    println(io, inner_lead*"+"*repeat("-", 70)*"+")
    
    if dm.selection !== nothing
        println(io, inner_lead*" ● Selection: Set")
        Base.show(io, dm.selection, vcat(level_code, 4))
    else
        println(io, inner_lead*" ○  Selection: Not Set")
    end
end