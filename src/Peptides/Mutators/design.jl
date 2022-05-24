"""
    DesignMutator(p_mut::AbstractFloat, grammar::LGrammar, selection::Opt{AbstractSelection}; [searchable_aminoacids::Dict{Char, Bool} = Peptides.available_aminoacids])

Return a [`DesignMutator`](@ref) instance. This `AbstractMutator` is a
_functor_, called with the following signature:

```
(design_mutator::DesignMutator)(pose::Pose)
```

The [`DesignMutator`](@ref) `AbstractMutator` loops through all [`Atom`](@ref)
instances in the given [`Pose`](@ref) and applies a mutation (using the
[`mutate!`](@ref ProtoSyn.Peptides.mutate!) method) if a random number
(`rand()`) is bellow a given probability of mutation `p_mut`. A different
mutation is attempted for each [`Atom`](@ref) instance looped over (therefore,
for a single attempt at a mutation per [`Residue`](@ref), unique [`Atom`](@ref)
names should be selected, as `an"CA"`, for example). If an `AbstractSelection`
`selection` is provided, only [`Atom`](@ref) instances marked as `true` in this
selection are considered for mutation. The applied mutation is sampled from the
given `searchable_aminoacids` dictionary (only for entries marked as `true`),
and uses the given [`LGrammar`](@ref) `grammar` to perform the mutation. If no
`searchable_aminoacids` dictionary is provided, will use the default
`Peptides.available_aminoacids` dictionary (contains all aminoacids set as
`true`). Note that the [`DesignMutator`](@ref) syncs any pending cartesian to
internal coordinate conversion (using the [`c2i!`](@ref ProtoSyn.c2i!) method).
Requests internal to cartesian coordinates conversion (using the
[`request_i2c!`](@ref ProtoSyn.request_i2c!) method). Does not [`sync!`](@ref)
the given [`Pose`](@ref) afterwards.

The [`DesignMutator`](@ref) `AbstractMutator` can also be optionally called
using the following signature, in which case only the provided list of
[`Atom`](@ref) instances will be considered for the application of this
`AbstractMutator`.

```
(design_mutator::DesignMutator)(pose::Pose, atoms::Vector{Atom})
```

# Fields
* `p_mut::AbtractFloat` - Compared against a `rand()` call, applies this Mutator to [`Atom`](@ref) instances where `rand() < p_mut`;
* `grammar::LGrammar` - The [`LGrammar`](@ref) instance with the aminoacid templates;
* `selection::Opt{AbstractSelection}` - If given, this Mutator will only loop over the selected [`Atom`](@ref) instances;
* `searchable_aminoacids::Dict{Char, Bool}` - The dictionary of available aminoacids, where each entry sets the availability of that aminoacid type to be the result of the mutation.

# See also
[`RotamerMutator`](@ref)

# Examples
```jldoctest
julia> dm = ProtoSyn.Peptides.Mutators.DesignMutator(1.0, res_lib, an"CA")
⚯  Design Mutator:
+----------------------------------------------------------------------+
| Index | Field                       | Value                          |
+----------------------------------------------------------------------+
| 1     | grammar                     | Set ✓                          |
| 2     | p_mut                       | 1.000                          |
+----------------------------------------------------------------------+
 ● Searchable aminoacids: M 🗹 | K 🗹 | P 🗹 | N 🗹 | I 🗹 | H 🗹 | E 🗹 | W 🗹 | S 🗹 | T 🗹 | C 🗹 | D 🗹 | A 🗹 | L 🗹 | Y 🗹 | V 🗹 | R 🗹 | G 🗹 | Q 🗹 | F 🗹 
 ● Selection: Set
 └── FieldSelection › Atom.name = CA
```
"""
mutable struct DesignMutator <: AbstractMutator
    p_mut::AbstractFloat
    grammar::LGrammar
    selection::Opt{AbstractSelection}
    searchable_aminoacids::Dict{Char, Bool}
end

DesignMutator(p_mut::AbstractFloat, grammar::LGrammar, selection::Opt{AbstractSelection}; searchable_aminoacids::Dict{Char, Bool} = deepcopy(Peptides.available_aminoacids)) = begin
    return DesignMutator(p_mut, grammar, selection, searchable_aminoacids)
end

function (design_mutator::DesignMutator)(pose::Pose)
    if design_mutator.selection === nothing
        atoms = Vector{Atom}(collect(eachatom(pose.graph)))
    else
        sele  = design_mutator.selection
        atoms = ProtoSyn.promote(sele, Atom)(pose, gather = true)
    end

    design_mutator(pose, atoms)
end

function (design_mutator::DesignMutator)(pose::Pose, atoms::Vector{Atom})
    searchable_aas = [a for (a, v) in design_mutator.searchable_aminoacids if v]

    # DesignMutator requires updated internal coordinates
    ProtoSyn.c2i!(pose.state, pose.graph) # Checks pose.state.c2i flag inside

    for atom in atoms
        if rand() < design_mutator.p_mut
            residue = atom.container

            # 1) Get different aminoacid
            cr_name = ProtoSyn.three_2_one[residue.name]
            nr_name = cr_name
            while nr_name == cr_name
                nr_name = rand(searchable_aas)
            end

            # 2) Perform mutation (already requests i2c)
            derivation = [string(nr_name)]
            Peptides.mutate!(pose, residue, design_mutator.grammar, derivation)
        end
    end
end

function Base.show(io::IO, dm::DesignMutator, level_code::Opt{LevelCode} = nothing)
    level_code = level_code === nothing ? LevelCode() : level_code
    lead       = ProtoSyn.get_lead(level_code)
    inner_lead = ProtoSyn.get_inner_lead(level_code)

    println(io, lead*"⚯  Design Mutator:")
    println(io, inner_lead*"+"*repeat("-", 70)*"+")
    @printf(io, "%s| %-5s | %-27s | %-30s |\n", inner_lead, "Index", "Field", "Value")
    println(io, inner_lead*"+"*repeat("-", 70)*"+")
    @printf(io, "%s| %-5d | %-27s | %-30s |\n", inner_lead, 1, "grammar", "Set ✓")
    @printf(io, "%s| %-5d | %-27s | %-30.3f |\n", inner_lead, 2, "p_mut", dm.p_mut)
    println(io, inner_lead*"+"*repeat("-", 70)*"+")

    s = ""
    for (aminoacid, searchable) = dm.searchable_aminoacids
        s *= searchable ? "$aminoacid 🗹 | " : "$aminoacid 🗷 | "
    end
    println(io, inner_lead*@sprintf " ● Searchable aminoacids: %s" s[1:(end-2)])
    
    if dm.selection !== nothing
        println(io, inner_lead*" ● Selection: Set")
        Base.show(io, dm.selection, vcat(level_code, 4))
    else
        println(io, inner_lead*" ○  Selection: Not Set")
    end
end