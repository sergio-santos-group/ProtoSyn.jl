export BondedToSelection
# Note: DistanceSelection is a BRANCH selection.

"""
    DistanceSelection(distance::Number, sele::S) where {S <: AbstractSelection}

A [`DistanceSelection`](@ref) takes an input selection `sele` and outputs a
[`Mask`](@ref) of [`Atom`](@ref) instances within the given `distance` (in
Ansgtrom Å) of the selected atoms from `sele`.

# State mode
The state mode of [`DistanceSelection`](@ref) `M` is forced to be Stateful

# Selection type

The selection type of [`DistanceSelection`](@ref) `T` is forced to be [`Atom`](@ref).

# Examples
```jldoctest
julia> sele = DistanceSelection(2.0, rn"ALA")
DistanceSelection ❯ Within 2.0 Å (Atom)
 └── FieldSelection › Residue.name = ALA

julia> 2.0:rn"ALA"
DistanceSelection ❯ Within 2.0 Å (Atom)
 └── FieldSelection › Residue.name = ALA
```
"""
mutable struct BondedToSelection{M, T} <: AbstractSelection
    sele::AbstractSelection

    BondedToSelection(sele::S) where {S <: AbstractSelection} = begin
        new{Stateless, Atom}(sele)
    end
end

state_mode_type(::BondedToSelection{M, T}) where {M, T} = M
selection_type(::BondedToSelection{M, T}) where {M, T} = T


# --- Select -------------------------------------------------------------------
function select(sele::BondedToSelection{Stateless, T}, container::AbstractContainer) where {T <: AbstractContainer}
    inner_mask = ProtoSyn.promote(sele.sele, Atom)(container)
    mask       = Mask{Atom}(falses(count_atoms(container)))
    count(inner_mask.content) === 0 && return mask
    atoms      = collect(eachatom(container))
    for (atom_i_index, atom_selected) in enumerate(inner_mask)
        !atom_selected && continue
        for bond in atoms[atom_i_index].bonds
            j = findfirst(x -> x === bond, atoms)
            mask[j] = true
        end
    end

    return mask
end

# --- Show ---------------------------------------------------------------------
function Base.show(io::IO, ds::BondedToSelection{M, T}, level_code::Opt{LevelCode} = nothing) where {M, T}
    lead = ProtoSyn.get_lead(level_code)
    if level_code === nothing
        level_code = LevelCode()
    end
    println(io, lead*"BondedToSelection ❯ ($T)")
    Base.show(io, ds.sele, vcat(level_code, 4))
end