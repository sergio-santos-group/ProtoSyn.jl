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
instances are looped over. Requests internal to cartesian coordinates conversion
(using [`request_i2c!`](@ref ProtoSyn.request_i2c!) method). Does not
[`sync!`](@ref) the given [`Pose`](@ref).

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
  Dihedral:
+--------------------------------------------------------------------+
| Index | Field                     | Value                          |
+--------------------------------------------------------------------+
| 1     | angle_sampler             | Function randn                 |
| 2     | p_mut                     | 1.0000                         |
| 3     | step_size                 | 1.0000                         |
| 4     | selection                 | Not set                        |
+--------------------------------------------------------------------+

julia> ProtoSyn.Mutators.DihedralMutator(randn, 0.05, 1.0, an"CA\$|C\$"r)
Dihedral:
+--------------------------------------------------------------------+
| Index | Field                     | Value                          |
+--------------------------------------------------------------------+
| 1     | angle_sampler             | Function randn                 |
| 2     | p_mut                     | 0.0500                         |
| 3     | step_size                 | 1.0000                         |
| 4     | selection                 | Set: FieldSelection            |
+--------------------------------------------------------------------+
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
        atoms = dihedral_mutator.selection(pose, gather = true)
    end
    
    dihedral_mutator(pose, atoms)
end

function (dihedral_mutator::DihedralMutator)(pose::Pose, atoms::Vector{Atom})
    for atom in atoms
        if rand() < dihedral_mutator.p_mut
            ∠ = dihedral_mutator.angle_sampler() * dihedral_mutator.step_size
            pose.state[atom].Δϕ += ∠
            ProtoSyn.request_i2c!(pose.state, all = true)
        end
    end
end

function Base.show(io::IO, dm::DihedralMutator)
    println("  Dihedral:")
    println(io, "+"*repeat("-", 68)*"+")
    @printf(io, "| %-5s | %-25s | %-30s |\n", "Index", "Field", "Value")
    println(io, "+"*repeat("-", 68)*"+")
    @printf(io, "| %-5d | %-25s | %-30s |\n", 1, "angle_sampler", "Function $(dm.angle_sampler)")
    @printf(io, "| %-5d | %-25s | %-30.4f |\n", 2, "p_mut", dm.p_mut)
    @printf(io, "| %-5d | %-25s | %-30.4f |\n", 3, "step_size", dm.step_size)
    if dm.selection === nothing
        @printf(io, "| %-5d | %-25s | %-30s |\n", 4, "selection", "Not set")
    else
        @printf(io, "| %-5d | %-25s | %-30s |\n", 4, "selection", "Set: $(typeof(dm.selection).name.name)")
    end
    println(io, "+"*repeat("-", 68)*"+")
end