using ProtoSyn.Peptides: Dihedral
using ProtoSyn: TrueSelection, getdihedral
using Printf

"""
    RotamerMutator(rotamer_library::Dict{String, ProtoSyn.Peptides.BBD_RotamerLibrary}, p_mut::AbstractFloat, n_first::Int, selection::Opt{AbstractSelection}

Return a [`RotamerMutator`](@ref) instance. This `AbstractMutator` is a
_functor_, called with the following signature:

```
(rotamer_mutator::RotamerMutator)(pose::Pose)
```

The [`RotamerMutator`](@ref) `AbstractMutator` loops through all [`Atom`](@ref)
instances in the given [`Pose`](@ref) and applies a [`Rotamer`](@ref ProtoSyn.Peptides.Rotamer)
conformation change if a random number (`rand()`) is bellow a given probability
of mutation `p_mut` (will skip any [`Residue`](@ref) with unnaccessible phi or
psi [`Dihedral`](@ref) angles, such as the first and last [`Residue`](@ref) of a
chain). A [`Rotamer`](@ref ProtoSyn.Peptides.Rotamer) conformation change is a concerted rotation of all
sidechain [`Atom`](@ref) instances in the [`Residue`](@ref) of the selected
[`Atom`](@ref) (therefore for a single attempt at [`Rotamer`](@ref ProtoSyn.Peptides.Rotamer) change,
unique [`Atom`](@ref) names should be selected, `an"CA"`, for example). If an
`AbstractSelection` `selection` is provided, only [`Atom`](@ref) instances
marked as `true` in this selection are considered for [`Rotamer`](@ref ProtoSyn.Peptides.Rotamer)
conformational change. The applied [`Rotamer`](@ref ProtoSyn.Peptides.Rotamer) is sampled from the
[`RotamerMutator`](@ref)`.rotamer_library`, based on the name of the
[`Residue`](@ref) and current phi and psi [`Dihedral`](@ref) angle values. The
`n_first` most likely [`Rotamer`](@ref ProtoSyn.Peptides.Rotamer) instances are taken into account during
this sampling step. Note that the [`RotamerMutator`](@ref) syncs any pending
cartesian to internal coordinate conversion (using the
[`c2i!`](@ref ProtoSyn.c2i!) method). Requests internal to cartesian coordinates
conversion (using [`request_i2c!`](@ref ProtoSyn.request_i2c!) method). Does not
[`sync!`](@ref) the given [`Pose`](@ref) afterwards.

The [`RotamerMutator`](@ref) `AbstractMutator` can also be optionally called
using the following signature, in which case only the provided list of
[`Atom`](@ref) instances will be considered for the application of this
`AbstractMutator`.

```
(rotamer_mutator::RotamerMutator)(pose::Pose, atoms::Vector{Atom})
```

# Fields
* `rotamer_library::Dict{String, ProtoSyn.Peptides.BBD_RotamerLibrary}` - A dictionary of [`BBD_RotamerLibrary`](@ref ProtoSyn.Peptides.BBD_RotamerLibrary) instances, for each aminoacid type;
* `p_mut::AbtractFloat` - Compared against a `rand()` call, applies this Mutator to [`Atom`](@ref) instances where `rand() < p_mut`;
* `n_first::Int` - Take only the N most likely [`Rotamer`](@ref ProtoSyn.Peptides.Rotamer) instances from the rotamer library;
* `selection::Opt{AbstractSelection}` - If given, this Mutator will only loop over the selected [`Atom`](@ref) instances;

# See also
[`DesignMutator`](@ref)

# Examples
```
julia> rm = ProtoSyn.Peptides.Mutators.RotamerMutator(rot_lib, 1.0, 10, an"CA" & !rn"PRO")
⚯  Rotamer Mutator:
+----------------------------------------------------------------------+
| Index | Field                       | Value                          |
+----------------------------------------------------------------------+
| 1     | rotamer_library             | Set ✓                          |
| 2     | p_mut                       | 1.000                          |
| 3     | n_first                     | 10                             |
+----------------------------------------------------------------------+
 ● Selection: Set
 └── BinarySelection ❯  & "and" (Atom)
      ├── FieldSelection › Atom.name = CA
      └── UnarySelection ❯ ! "not" (Residue)
           └── FieldSelection › Residue.name = PRO
```
"""
mutable struct RotamerMutator <: AbstractMutator
    rotamer_library::Dict{String, ProtoSyn.Peptides.BBD_RotamerLibrary}
    p_mut::AbstractFloat
    n_first::Int
    selection::Opt{AbstractSelection}
end


function (rotamer_mutator::RotamerMutator)(pose::Pose)
    if rotamer_mutator.selection === nothing
        atoms = Vector{Atom}(collect(eachatom(pose.graph)))
    else
        sele = rotamer_mutator.selection
        atoms = ProtoSyn.promote(sele, Atom)(pose, gather = true)
    end

    rotamer_mutator(pose, atoms)
end


function (rotamer_mutator::RotamerMutator)(pose::Pose, atoms::Vector{Atom})
    
    # RotamerMutator requires updated internal coordinates
    ProtoSyn.c2i!(pose.state, pose.graph) # Checks pose.state.c2i flag inside
    
    for atom in atoms
        if rand() < rotamer_mutator.p_mut

            # 1) Get the residue name
            residue = atom.container
            name    = residue.name

            # 2) Get the residue phi & psi
            phi_atom = Dihedral.phi(residue)
            phi_atom === nothing && continue
            phi = getdihedral(pose.state, phi_atom)
            psi_atom = Dihedral.psi(residue)
            psi_atom === nothing && continue
            psi = getdihedral(pose.state, psi_atom)

            # 3) Sample with correct name, phi & psi
            rotamer_stack = nothing
            try
                rotamer_stack = rotamer_mutator.rotamer_library[name][phi, psi]
            catch KeyError
                continue
            end
            rotamer = Peptides.sample(rotamer_stack, rotamer_mutator.n_first)

            # 4) Apply sampled rotamer
            Peptides.apply!(pose.state, rotamer, residue)
            ProtoSyn.request_i2c!(pose.state, all = true)
        end
    end
end

function Base.show(io::IO, rm::RotamerMutator, level_code::Opt{LevelCode} = nothing)
    level_code = level_code === nothing ? LevelCode() : level_code
    lead       = ProtoSyn.get_lead(level_code)
    inner_lead = ProtoSyn.get_inner_lead(level_code)

    println(io, lead*"⚯  Rotamer Mutator:")
    println(io, inner_lead*"+"*repeat("-", 70)*"+")
    @printf(io, "%s| %-5s | %-27s | %-30s |\n", inner_lead, "Index", "Field", "Value")
    println(io, inner_lead*"+"*repeat("-", 70)*"+")
    @printf(io, "%s| %-5d | %-27s | %-30s |\n", inner_lead, 1, "rotamer_library", "Set ✓")
    @printf(io, "%s| %-5d | %-27s | %-30.3f |\n", inner_lead, 2, "p_mut", rm.p_mut)
    @printf(io, "%s| %-5d | %-27s | %-30d |\n", inner_lead, 3, "n_first", rm.n_first)
    println(io, inner_lead*"+"*repeat("-", 70)*"+")
    
    if rm.selection !== nothing
        println(io, inner_lead*" ● Selection: Set")
        Base.show(io, rm.selection, vcat(level_code, 4))
    else
        println(io, inner_lead*" ○  Selection: Not Set")
    end
end