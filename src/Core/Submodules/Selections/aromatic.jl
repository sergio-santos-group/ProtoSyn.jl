export AromaticSelection
# Note: AromaticSelection is a LEAF selection.

"""
    AromaticSelection([d::Int = 6])
    

An [`AromaticSelection`](@ref) selects [`Atom`](@ref) instances that belong to
an aromatic ring. Rings are identified by their bonds, using the
[`travel_bonds`](@ref) method: any atom who exists in the set of bonded atoms
within `d` bonds to itself forms a ring (d = 6, by default). 

# State mode

The state mode of [`SerialSelection`](@ref) `M` is forced to be `Stateless`.

# Selection type

The selection type of [`SerialSelection`](@ref) `T` is forced to be
[`Atom`](@ref).

!!! ukw "Note:"
    This selection does not have a short syntax version.

# Examples
```jldoctest
julia> AromaticSelection(6)
AromaticSelection (Rings with 6 or less bonds) › (Atom)

julia> AromaticSelection(7)
AromaticSelection (Rings with 7 or less bonds) › (Atom)
```
"""
mutable struct AromaticSelection{M, T} <: AbstractSelection

    d::Int

    AromaticSelection(d::Int = 6) = begin
        new{ProtoSyn.Stateless, Atom}(d)
    end
end

state_mode_type(::AromaticSelection{M, T}) where {M, T} = M
selection_type(::AromaticSelection{M, T}) where {M, T} = T

# --- Select -------------------------------------------------------------------
function select(sele::AromaticSelection, container::AbstractContainer)

    n_items = count_atoms(container)
    mask    = Mask{Atom}(n_items)

    for (index, atom) in enumerate(eachatom(container))
        range = ProtoSyn.travel_bonds(atom, sele.d)
        if atom in range[2:end]
            mask[index] = true
        end
    end
    return mask
end

# --- Show ---------------------------------------------------------------------
function Base.show(io::IO, as::AromaticSelection{M, T}, level_code::Opt{LevelCode} = nothing) where {M, T}
    lead = ProtoSyn.get_lead(level_code)
    if level_code === nothing
        level_code = LevelCode()
    end
    println(io, lead*"AromaticSelection (Rings with $(as.d) or less bonds) › ($(selection_type(as)))")
end