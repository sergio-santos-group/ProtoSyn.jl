export ChargeSelection
# Note: ChargeSelection is a LEAF selection.

"""
    ChargeSelection(charge::F, op::Function = ===) where {F <: AbstractFloat}

A [`ChargeSelection`](@ref) takes an input `charge` value and an `op` (as a
`Function`, is `===`, by default) and outputs a `Mask` (of type
[`Atom`](@ref)) containing all [`Atom`](@ref) instances whose `AtomState` charge
matches the given `charge` (when compared with `op`). As an example, `op` can be
`>`, `<=`, `==`, etc.

# State mode

The state mode of [`ChargeSelection`](@ref) `M` is forced to be `Stateful`.

# Selection type

The selection type of [`ChargeSelection`](@ref) is forced to be [`Atom`](@ref).

!!! ukw "Note:"
    This selection does not have a short syntax version.

# Examples
```jldoctest
julia> ChargeSelection(1.0)
ChargeSelection (Atoms with charge === 1.0) › (Atom)

julia> ChargeSelection(-1.0, >)
ChargeSelection (Atoms with charge > -1.0) › (Atom)
```
"""
mutable struct ChargeSelection{M, T} <: AbstractSelection

    charge::Float64
    op::Function

    ChargeSelection(charge::Float64, op::Function = ===) = begin
        new{ProtoSyn.Stateful, Atom}(charge, op)
    end
end

state_mode_type(::ChargeSelection{M, T}) where {M, T} = M
selection_type(::ChargeSelection{M, T}) where {M, T} = T

# --- Select -------------------------------------------------------------------
function select(sele::ChargeSelection{Stateful, T}, container::AbstractContainer) where {T <: AbstractContainer}

    n_items = counter(T)(container)
    mask    = Mask{T}(n_items)

    return function (state::State)
        for (index, atom) in enumerate(eachatom(container))
            if sele.op(state[atom].δ, sele.charge)
                mask[index] = true
            end
        end

        return mask
    end
end

# --- Show ---------------------------------------------------------------------
function Base.show(io::IO, cs::ChargeSelection{M, T}, level_code::Opt{LevelCode} = nothing) where {M, T}
    lead = ProtoSyn.get_lead(level_code)
    if level_code === nothing
        level_code = LevelCode()
    end
    println(io, lead*"ChargeSelection (Atoms with charge $(cs.op) $(cs.charge)) › ($(selection_type(cs)))")
end