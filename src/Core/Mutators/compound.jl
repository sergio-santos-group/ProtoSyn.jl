"""
    CompoundMutator(mutators::Vector{AbstractMutator}, selection::Union{AbstractSelection, Nothing})

Return a [`CompoundMutator`](@ref) instance. This `AbstractMutator` is a
_functor_, called with the following signature:

```
(compound_mutator::CompoundMutator)(pose::Pose)
```

A [`CompoundMutator`](@ref) loops over a list of one or more `AbstractMutator`
instances, applying them in the provided sequence to a [`Pose`](@ref) instance.
If an `AbstractSelection` `selection` is provided, the same selection will be
used when applying each `AbstractMutator`, merged with any possible inner
selection specific for each mutator (using the & "and" operator). For example, a
[`DihedralMutator`](@ref) may have an `an"CA"` selection attached. When called
from a [`CompoundMutator`](@ref) with a `rid"1:10"` selection attached, only the
`CA` [`Atom`](@ref) instances from the region between [`Residue`](@ref) 1 to
[`Residue`](@ref) 10 will be considered for application of the
[`DihedralMutator`](@ref). Note that, when applying each mutator, any required
[`sync!`](@ref) call and request for internal to cartesian coordinate
conversion (or vice-versa) are handled by each `AbstractMutator` independently.

The [`CompoundMutator`](@ref) `AbstractMutator` can also be optionally called
using the following signature, in which case only the provided list of
[`Atom`](@ref) instances will be considered for the application of this
`AbstractMutator` (i.e.: [`CompoundMutator`](@ref) of a
[`CompoundMutator`](@ref)).

```
(compound_mutator::CompoundMutator)(pose::Pose, atoms::Vector{Atom})
```

!!! ukw "Note:"
    This `AbstractMutator` requires that the mutators provided in the list have
    a `(compound_mutator::CompoundMutator)(pose::Pose, atoms::Vector{Atom})`
    signature. Check each `AbstractMutator` documentation entry and
    [Creating custom mutators](@ref) section for details. 

# Fields
* `mutators::Vector{AbstractMutator}` - The list of `AbstractMutator` instances this mutator will loop over and apply sequentially;
* `selection::Union{AbstractSelection, Nothing}` - If provided, all inner `AbstractMutator` instances will be applied to each individual selection (if present) merged with this common selection.

# See also
[`CompoundDriver`](@ref)

# Examples
```jldoctest
julia> m = ProtoSyn.Mutators.CompoundMutator([rrbm, trbm], sele)
  Compound Mutator:
   ↘ Mutators:

  Rotation Rigid Body Mutator:
+----------------------------------------------------------------------+
| Index | Field                       | Value                          |
+----------------------------------------------------------------------+
| 1     | axis_sampler                | Function rand_vector_in_sphere |
| 2     | angle_sampler               | Function randn                 |
| 3     | pivot_sampler               | Function center_of_mass        |
| 4     | step_size                   | 0.4000                         |
| 5     | selection                   | Not set                        |
+----------------------------------------------------------------------+

  Translation Rigid Body Mutator:
+----------------------------------------------------------------------+
| Index | Field                       | Value                          |
+----------------------------------------------------------------------+
| 1     | translation_vector_sampler  | Function rand_vector_in_sphere |
| 2     | step_size                   | 0.4000                         |
| 3     | selection                   | Not set                        |
+----------------------------------------------------------------------+

   ↘ Selection:
RandomSelectionFromList ❯ (Residue)
 ├── RangeSelection › Residue.id between 1 and 22
 └── RangeSelection › Residue.id between 25 and 46
```

"""
mutable struct CompoundMutator <: AbstractMutator
    mutators::Vector{AbstractMutator}
    selection::Union{AbstractSelection, Nothing}
end

function (compound_mutator::CompoundMutator)(pose::Pose)

    if compound_mutator.selection !== nothing
        sele  = ProtoSyn.promote(compound_mutator.selection, Atom)
    else
        sele = TrueSelection{Atom}()
    end
    atoms = sele(pose, gather = true)

    # Note: From one mutator to the next, certain changes may have taken place.
    # If these changes are applied to the internal coordinates only, no sync!
    # step is performed. Example: 2 DihedralMutators in a row. However, each 
    # Mutator should internally verify that the required set of coordinates has
    # been correctly updated. Example: CrankshaftMutators check i2c! since they
    # require the updated cartesian coordinates.
    for mutator in compound_mutator.mutators
        mutator(pose, atoms)
    end
end


function (compound_mutator::CompoundMutator)(pose::Pose, atoms::Vector{Atom})

    if compound_mutator.selection !== nothing
        sele  = ProtoSyn.promote(compound_mutator.selection, Atom)
    else
        sele = TrueSelection{Atom}()
    end
    _atoms = sele(pose, gather = true)
    atoms = intersect(_atoms, atoms)

    # Note: From one mutator to the next, certain changes may have taken place.
    # If these changes are applied to the internal coordinates only, no sync!
    # step is performed. Example: 2 DihedralMutators in a row. However, each 
    # Mutator should internally verify that the required set of coordinates has
    # been correctly updated. Example: CrankshaftMutators check i2c! since they
    # require the updated cartesian coordinates.
    for mutator in compound_mutator.mutators
        mutator(pose, atoms)
    end
end


function Base.show(io::IO, cm::CompoundMutator)
    println(io, "  Compound Mutator:\n   ↘ Mutators:\n")
    for mutator in cm.mutators
        println(io, mutator)
    end
    println(io, "   ↘ Selection:")
    println(io, cm.selection)

end