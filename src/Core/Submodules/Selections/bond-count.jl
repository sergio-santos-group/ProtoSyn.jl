export BondCountSelection
# Note: BondCountSelection is a LEAF selection.

"""
    BondCountSelection{T}(n::Int, op::Function = ===)

A [`BondCountSelection`](@ref) takes an input `n` value and an `op` (as a
`Function`, is `===`, by default) and outputs a `Mask` (of type
[`Atom`](@ref)) containing all [`Atom`](@ref) instances whose number of bonded
[`Atom`](@ref) instances (in `atom.bonds`) matches the given `n` (when compared
with `op`). As an example, `op` can be `>`, `<=`, `==`, etc.

# State mode

The state mode of [`BondCountSelection`](@ref) `M` is forced to be `Stateless`.

# Selection type

The selection type of [`BondCountSelection`](@ref) is forced to be [`Atom`](@ref).

!!! ukw "Note:"
    This selection does not have a short syntax version.

# Examples
```jldoctest
julia> BondCountSelection(1)
BondCountSelection (Atoms with N bonds === 1) › (Atom)

julia> BondCountSelection(3, >)
BondCountSelection (Atoms with N bonds > 3) › (Atom)
```
"""
mutable struct BondCountSelection{M, T} <: AbstractSelection

    n::Int
    op::Function

    BondCountSelection(n::Int, op::Function = ===) = begin
        new{ProtoSyn.Stateless, Atom}(n, op)
    end
end

state_mode_type(::BondCountSelection{M, T}) where {M, T} = M
selection_type(::BondCountSelection{M, T}) where {M, T} = T

# --- Select -------------------------------------------------------------------
function select(sele::BondCountSelection{Stateless, T}, container::AbstractContainer) where {T <: AbstractContainer}

    n_items = counter(T)(container)
    mask    = Mask{T}(n_items)

    for (index, atom) in enumerate(eachatom(container))
        if sele.op(length(atom.bonds), sele.n)
            mask[index] = true
        end
    end

    return mask
end

# --- Show ---------------------------------------------------------------------
function Base.show(io::IO, bcs::BondCountSelection{M, T}, level_code::Opt{LevelCode} = nothing) where {M, T}
    lead = ProtoSyn.get_lead(level_code)
    if level_code === nothing
        level_code = LevelCode()
    end
    println(io, lead*"BondCountSelection (Atoms with N bonds $(bcs.op) $(bcs.n)) › ($(selection_type(bcs)))")
end