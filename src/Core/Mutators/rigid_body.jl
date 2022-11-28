# ------------------------------------------------------------------------------
# Translation Rigid Body

"""
    TranslationRigidBodyMutator(translation_vector_sampler::Function, step_size::AbstractFloat, selection::Opt{AbstractSelection})

Return a new [`TranslationRigidBodyMutator`](@ref) instance. This
`AbstractMutator` is a _functor_, called with the following signature:

```
(rigid_body_mutator::TranslationRigidBodyMutator)(pose::Pose)
```

The [`TranslationRigidBodyMutator`](@ref) `AbstractMutator` applies a
translation to all [`Atom`](@ref) instances in a given [`Pose`](@ref) based on a
given axis. This axis is sampled from `translation_vector_sampler` (who receives
no input arguments and should return a `Vector{Float}` with size `3`, the `X`,
`Y` and `Z` dimensions). Although not necessary, this `Vector{Float}` should
have norm `1.0`. The translation vector is then multiplied by `step_size`. If an
`AbstractSelection` `selection` is provided, only the selected [`Atom`](@ref)
instances are translated. Note that the [`TranslationRigidBodyMutator`](@ref)
syncs any pending internal to cartesian coordinate conversion (using the
[`i2c!`](@ref ProtoSyn.i2c!) method). Requests cartesian to internal coordinates
conversion (using [`request_c2i!`](@ref ProtoSyn.request_c2i!) method). Does not
[`sync!`](@ref) the given [`Pose`](@ref) afterwards.

The [`TranslationRigidBodyMutator`](@ref) `AbstractMutator` can also be
optionally called using the following signature, in which case only the provided
list of [`Atom`](@ref) instances will be considered for the application of this
`AbstractMutator`.

```
(rigid_body_mutator::TranslationRigidBodyMutator)(pose::Pose, atoms::Vector{Atom})
```

# Fields

* `translation_vector_sampler::Function` - Should return a `Vector{Float}` axis (X, Y and Z dimensions). Is called with no input arguments;
* `step_size::AbstractFloat` - Multiplies the sampled axis by this value;
* `selection::Opt{AbstractSelection}` - If given, this Mutator will only translate the selected [`Atom`](@ref) instances;

# See also
[`TranslationRigidBodyMutator`](@ref)

# Examples
```jldoctest
julia> m = ProtoSyn.Mutators.TranslationRigidBodyMutator(ProtoSyn.rand_vector_in_sphere, 1.0, rn"CBZ")
⚯  Translation Rigid Body Mutator:
+----------------------------------------------------------------------+
| Index | Field                       | Value                          |
+----------------------------------------------------------------------+
| 1     | translation_vector_sampler  | Function rand_vector_in_sphere |
| 2     | step_size                   | 1.0000                         |
+----------------------------------------------------------------------+
 ● Selection: Set
 └── FieldSelection › Residue.name = CBZ
```
"""
mutable struct TranslationRigidBodyMutator <: AbstractMutator
    translation_vector_sampler::Function # Should return a Vector{Float}
    step_size::AbstractFloat
    selection::Opt{AbstractSelection}
end


function (rigid_body_mutator::TranslationRigidBodyMutator)(pose::Pose)
    if rigid_body_mutator.selection === nothing
        atoms = Vector{Atom}(collect(eachatom(pose.graph)))
    else
        sele  = rigid_body_mutator.selection
        atoms = ProtoSyn.promote(sele, Atom)(pose, gather = true)
    end
    
    rigid_body_mutator(pose, atoms)
end


function (rigid_body_mutator::TranslationRigidBodyMutator)(pose::Pose, atoms::Vector{Atom})

    # TranslationRigidBodyMutator requires updated cartesian coordinates
    ProtoSyn.i2c!(pose.state, pose.graph) # Checks pose.state.i2c flag inside

    translation_vector   = rigid_body_mutator.translation_vector_sampler()
    translation_vector .*= rigid_body_mutator.step_size
    for atom in atoms
        # Note: Since the :t setproperty! is being intercepted, we can't do
        # pose.state[atom].t .+= translation_vector
        pose.state[atom].t = translation_vector .+ pose.state[atom].t
    end
    ProtoSyn.request_c2i!(pose.state)
end

function Base.show(io::IO, trbm::TranslationRigidBodyMutator, level_code::Opt{LevelCode} = nothing)
    lead = ProtoSyn.get_lead(level_code)
    if level_code === nothing
        level_code = LevelCode()
        inner_lead = lead
    else
        inner_level_code = copy(level_code)
        inner_level_code[end] = level_code.conv_table[level_code.levels[end]]
        inner_lead = ProtoSyn.get_lead(inner_level_code)
    end


    println(io, lead*"⚯  Translation Rigid Body Mutator:")
    println(io, inner_lead*"+"*repeat("-", 80)*"+")
    @printf(io, "%s| %-5s | %-27s | %-40s |\n", inner_lead, "Index", "Field", "Value")
    println(io, inner_lead*"+"*repeat("-", 80)*"+")
    @printf(io, "%s| %-5d | %-27s | %-40s |\n", inner_lead, 1, "translation_vector_sampler", "Function $(trbm.translation_vector_sampler)")
    @printf(io, "%s| %-5d | %-27s | %-40.4f |\n", inner_lead, 2, "step_size", trbm.step_size)
    println(io, inner_lead*"+"*repeat("-", 80)*"+")
    
    if trbm.selection !== nothing
        println(io, inner_lead*" ● Selection: Set")
        Base.show(io, trbm.selection, vcat(level_code, 4))
    else
        println(io, inner_lead*" ○  Selection: Not Set")
    end
end

# ------------------------------------------------------------------------------
# Rotation Rigid Body

"""
    RotationRigidBodyMutator(axis_sampler::Function, angle_sampler::Function, pivot_sampler::Function, step_size::AbstractFloat, selection::Opt{AbstractSelection})

Return a new [`RotationRigidBodyMutator`](@ref) instance. This `AbstractMutator`
is a _functor_, called with the following signature:

```
(rigid_body_mutator::RotationRigidBodyMutator)(pose::Pose)
```

The [`RotationRigidBodyMutator`](@ref) `AbstractMutator` applies a rotation to
all [`Atom`](@ref) instances in a given [`Pose`](@ref) based on a given angle
and axis, and centered around a pivot position (See
[`rotation_matrix_from_axis_angle`](@ref ProtoSyn.rotation_matrix_from_axis_angle)).
The applied rotation's angle is sampled by calling `angle_sampler` (receives no
input arguments and should return a `Float` angle value in radians) and
multiplied by the `step_size` value. The rotation's axis is sampled from the
`axis_sampler` (receives no input arguments and should return a `Vector{Float}`
with size `3`, the `X`, `Y` and `Z` dimensions). Finally, this rotation is
centered in a pivot position, sampled from `pivot_sampler` (receives the
[`Pose`](@ref) and a list of selected [`Atom`](@ref) indexes as the input
arguments and should return a `Vector{Float}` with size `3`, the `X`, `Y` and
`Z` dimensions). If an `AbstractSelection` `selection` is provided, only the
selected [`Atom`](@ref) instances are rotated and used as input for
`pivot_sampler`. Note that the [`RotationRigidBodyMutator`](@ref) syncs any
pending internal to cartesian coordinate conversion (using the
[`i2c!`](@ref ProtoSyn.i2c!) method). Requests cartesian to internal coordinates
conversion (using [`request_c2i!`](@ref ProtoSyn.request_c2i!) method). Does not
[`sync!`](@ref) the given [`Pose`](@ref) afterwards.

The [`RotationRigidBodyMutator`](@ref) `AbstractMutator` can also be
optionally called using the following signature, in which case only the provided
list of [`Atom`](@ref) instances will be considered for the application of this
`AbstractMutator`.

```
(rigid_body_mutator::RotationRigidBodyMutator)(pose::Pose, atoms::Vector{Atom})
```

# Fields

* `axis_sampler::Function` - Should return a `Vector{Float}` axis (X, Y and Z dimensions). Is called with no input arguments;
* `angle_sampler::Function` - Should return a `Float` angle value (in radians). Is called with no input arguments;
* `pivot_sampler::Function` - Should return a `Vector{Float}` position (X, Y and Z dimensions). Is called with 2 input arguments: a [`Pose`](@ref) instance and a `Vector{Int}` with the indexes of the selected [`Atom`](@ref) instances;
* `step_size::AbstractFloat` - Multiplies the sampled angle by this value;
* `selection::Opt{AbstractSelection}` - If given, this Mutator will only be applied to the selected [`Atom`](@ref) instances.

# See also
[`TranslationRigidBodyMutator`](@ref)

# Examples
```jldoctest
julia> m = ProtoSyn.Mutators.RotationRigidBodyMutator(ProtoSyn.rand_vector_in_sphere, randn, ProtoSyn.center_of_mass, 1.0, rn"CBZ")
⚯  Rotation Rigid Body Mutator:
+----------------------------------------------------------------------+
| Index | Field                       | Value                          |
+----------------------------------------------------------------------+
| 1     | axis_sampler                | Function rand_vector_in_sphere |
| 2     | angle_sampler               | Function randn                 |
| 3     | pivot_sampler               | Function center_of_mass        |
| 4     | step_size                   | 1.0000                         |
+----------------------------------------------------------------------+
 ● Selection: Set
 └── FieldSelection › Residue.name = CBZ
```
"""
mutable struct RotationRigidBodyMutator <: AbstractMutator
    axis_sampler::Function # Should return a Vector{Float}
    angle_sampler::Function # Should return a Float
    pivot_sampler::Function # Receives a pose and a vector of selected
    # atom indexes and should return a Vector{Float}.
    step_size::AbstractFloat
    selection::Opt{AbstractSelection}
end


function (rigid_body_mutator::RotationRigidBodyMutator)(pose::Pose)

    # RotationRigidBodyMutator requires updated cartesian coordinates
    sync!(pose)
    # ProtoSyn.i2c!(pose.state, pose.graph) # Checks pose.state.i2c flag inside

    axis  = rigid_body_mutator.axis_sampler()
    angle = rigid_body_mutator.angle_sampler() * rigid_body_mutator.step_size
    rmat  = ProtoSyn.rotation_matrix_from_axis_angle(axis, angle)
    selec = rigid_body_mutator.selection === nothing ? TrueSelection{Atom}() : rigid_body_mutator.selection
    mask  = ProtoSyn.promote(selec, Atom)(pose)
    idxs  = findall(mask.content)
    pivot = rigid_body_mutator.pivot_sampler(pose, idxs)
    pose.state.x[:, idxs] = (rmat * (pose.state.x[:, idxs] .- pivot)) .+ pivot
    ProtoSyn.request_c2i!(pose.state)
end

function (rigid_body_mutator::RotationRigidBodyMutator)(pose::Pose, atoms::Vector{Atom})

    # RotationRigidBodyMutator requires updated cartesian coordinates
    sync!(pose)
    # ProtoSyn.i2c!(pose.state, pose.graph) # Checks pose.state.i2c flag inside

    axis  = rigid_body_mutator.axis_sampler()
    angle = rigid_body_mutator.angle_sampler() * rigid_body_mutator.step_size
    rmat  = ProtoSyn.rotation_matrix_from_axis_angle(axis, angle)
    idxs = ProtoSyn.ids(atoms)
    pivot = rigid_body_mutator.pivot_sampler(pose, idxs)
    pose.state.x[:, idxs] = (rmat * (pose.state.x[:, idxs] .- pivot)) .+ pivot
    ProtoSyn.request_c2i!(pose.state)
end

function Base.show(io::IO, rrbm::RotationRigidBodyMutator, level_code::Opt{LevelCode} = nothing)
    lead = ProtoSyn.get_lead(level_code)

    if level_code === nothing
        level_code = LevelCode()
        inner_lead = lead
    else
        inner_level_code = copy(level_code)
        inner_level_code[end] = level_code.conv_table[level_code.levels[end]]
        inner_lead = ProtoSyn.get_lead(inner_level_code)
    end

    println(io, lead*"⚯  Rotation Rigid Body Mutator:")
    println(io, inner_lead*"+"*repeat("-", 80)*"+")
    @printf(io, "%s| %-5s | %-27s | %-40s |\n", inner_lead, "Index", "Field", "Value")
    println(io, inner_lead*"+"*repeat("-", 80)*"+")
    @printf(io, "%s| %-5d | %-27s | %-40s |\n", inner_lead, 1, "axis_sampler", "Function $(rrbm.axis_sampler)")
    @printf(io, "%s| %-5d | %-27s | %-40s |\n", inner_lead, 2, "angle_sampler", "Function $(rrbm.angle_sampler)")
    @printf(io, "%s| %-5d | %-27s | %-40s |\n", inner_lead, 3, "pivot_sampler", "Function $(rrbm.pivot_sampler)")
    @printf(io, "%s| %-5d | %-27s | %-40.4f |\n", inner_lead, 4, "step_size", rrbm.step_size)
    println(io, inner_lead*"+"*repeat("-", 80)*"+")
    
    if rrbm.selection !== nothing
        println(io, inner_lead*" ● Selection: Set")
        Base.show(io, rrbm.selection, vcat(level_code, 4))
    else
        println(io, inner_lead*" ○  Selection: Not Set")
    end
end