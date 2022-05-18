export AromaticSelection
# Note: FieldSelection is a LEAF selection.

"""
# TODO: Documentation
"""
mutable struct AromaticSelection{M, T} <: AbstractSelection
    AromaticSelection() = begin
        new{ProtoSyn.Stateless, Atom}()
    end
end

state_mode_type(::AromaticSelection{M, T}) where {M, T} = M
selection_type(::AromaticSelection{M, T}) where {M, T} = T

# --- Select -------------------------------------------------------------------
function select(sele::AromaticSelection, container::AbstractContainer)

    n_items = count_atoms(container)
    mask = Mask{Atom}(n_items)

    for (index, atom) in enumerate(eachatom(container))
        range = ProtoSyn.travel_bonds(atom, 6)
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
    println(io, lead*"AromaticSelection › ($(selection_type(as)))")
end