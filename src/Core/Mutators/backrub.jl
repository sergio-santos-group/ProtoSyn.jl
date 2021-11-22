"""
    BackrubMutator(translation_vector_sampler::Function, p_mut::AbstractFloat, step_size::AbstractFloat, selection::Opt{AbstractSelection})

Return a [`BackrubMutator`](@ref) instance. This `AbstractMutator` is a
_functor_, called with the following signature:

```
(backrub_mutator::BackrubMutator)(pose::Pose)
```

The [`BackrubMutator`](@ref) `AbstractMutator` loops through all the
[`Atom`](@ref) instances in a given [`Pose`](@ref) and applies a translation
movement if a random number (`rand()`) is bellow a given probability of mutation
`p_mut` (therefore a higher `p_mut` value applies a larger number of translation
movements per call). The translation vector is sampled from the
`translation_vector_sampler`, multiplied by the given `step_size`. The
resulting value is then added to the selected [`Atom`](@ref)`.t`. Note that a
new translation vector is sampled for each selected [`Atom`](@ref) instance. If
an `AbstractSelection` `selection` is provided, only the selected [`Atom`](@ref)
instances are looped over. If the given `AbstractSelection` `selection` is not
of selection type [`Atom`](@ref), it will be promoted to this type (using
[`promote`](@ref ProtoSyn.promote) with default aggregator `any`). Note that the
[`BackrubMutator`](@ref) syncs any pending internal to cartesian coordinate
conversion (using the [`i2c!`](@ref ProtoSyn.i2c!) method). Requests cartesian
to internal coordinates conversion (using
[`request_c2i!`](@ref ProtoSyn.request_c2i!) method). Does not [`sync!`](@ref)
the given [`Pose`](@ref) afterwards.

The [`BackrubMutator`](@ref) `AbstractMutator` can also be optionally called
using the following signature, in which case only the provided list of
[`Atom`](@ref) instances will be considered for the application of this
`AbstractMutator`.

```
(backrub_mutator::BackrubMutator)(pose::Pose, atoms::Vector{Atom})
```

# Fields

* `translation_vector_sampler::Function` - Should return a `Vector{Float}` axis (X, Y and Z dimensions). Is called with no input arguments;
* `p_mut::AbtractFloat` - Compared against a `rand()` call, applies this Mutator to [`Atom`](@ref) instances where `rand() < p_mut`;
* `step_size::AbstractFloat` - Multiplies the sampled vector by this value;
* `selection::Opt{AbstractSelection}` - If given, this Mutator will only loop over the selected [`Atom`](@ref) instances.

# See also
[`TranslationRigidBodyMutator`](@ref)

# Examples
```jldoctest
⚯  Backrub Mutator:
+----------------------------------------------------------------------+
| Index | Field                       | Value                          |
+----------------------------------------------------------------------+
| 1     | translation_vector_sampler  | Function rand_vector_in_sphere |
| 2     | p_mut                       | 1.0000                         |
| 3     | step_size                   | 1.0000                         |
+----------------------------------------------------------------------+
 ○  Selection: Not Set
```
"""
mutable struct BackrubMutator <: AbstractMutator
    translation_vector_sampler::Function # Should return a Vector{Float}
    p_mut::AbstractFloat
    step_size::AbstractFloat
    selection::Opt{AbstractSelection}
end


function (backrub_mutator::BackrubMutator)(pose::Pose)
    if backrub_mutator.selection === nothing
        atoms = Vector{Atom}(collect(eachatom(pose.graph)))
    else
        if ProtoSyn.selection_type(backrub_mutator.selection) !== Atom
            sele  = ProtoSyn.promote(backrub_mutator.selection, Atom) # default aggregator is `any`
            atoms = sele(pose, gather = true)
        else
            atoms = backrub_mutator.selection(pose, gather = true)
        end
    end
    
    backrub_mutator(pose, atoms)
end


function (backrub_mutator::BackrubMutator)(pose::Pose, atoms::Vector{Atom})

    # BackrubMutator requires updated cartesian coordinates
    ProtoSyn.i2c!(pose.state, pose.graph) # Checks pose.state.i2c flag inside

    for atom in atoms
        if rand() < backrub_mutator.p_mut
            translation_vector   = backrub_mutator.translation_vector_sampler()
            translation_vector .*= backrub_mutator.step_size
            pose.state[atom].t   = translation_vector .+ pose.state[atom].t
        end
    end

    ProtoSyn.request_c2i!(pose.state)
end


function Base.show(io::IO, bm::BackrubMutator, level_code::Opt{LevelCode} = nothing)
    lead = ProtoSyn.get_lead(level_code)
    if level_code === nothing
        level_code = LevelCode()
        inner_lead = lead
    else
        inner_level_code = copy(level_code)
        inner_level_code[end] = level_code.conv_table[level_code.levels[end]]
        inner_lead = ProtoSyn.get_lead(inner_level_code)
    end


    println(io, lead*"⚯  Backrub Mutator:")
    println(io, inner_lead*"+"*repeat("-", 70)*"+")
    @printf(io, "%s| %-5s | %-27s | %-30s |\n", inner_lead, "Index", "Field", "Value")
    println(io, inner_lead*"+"*repeat("-", 70)*"+")
    @printf(io, "%s| %-5d | %-27s | %-30s |\n", inner_lead, 1, "translation_vector_sampler", "Function $(bm.translation_vector_sampler)")
    @printf(io, "%s| %-5d | %-27s | %-30.4f |\n", inner_lead, 2, "p_mut", bm.p_mut)
    @printf(io, "%s| %-5d | %-27s | %-30.4f |\n", inner_lead, 3, "step_size", bm.step_size)
    println(io, inner_lead*"+"*repeat("-", 70)*"+")
    
    if bm.selection !== nothing
        println(io, inner_lead*" ● Selection: Set")
        Base.show(io, bm.selection, vcat(level_code, 4))
    else
        println(io, inner_lead*" ○  Selection: Not Set")
    end
end