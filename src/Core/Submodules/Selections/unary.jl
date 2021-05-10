export UnarySelection
# Note: UnarySelection is a BRANCH selection.

"""
    UnarySelection{M}(op::Function, sele::AbstractSelection)

A [`UnarySelection`](@ref) applies an operation `op` to the given `AbstractSelection`
`sele`. Available operations with short syntax:
* `!` - Logical Not. Negates the values selected by `sele`.

# State mode

The selection type of [`UnarySelection`](@ref) can be either `Stateless` or
`Stateful`. When using the short syntax, it will automatically be set to the
`StateMode` of the provided `sele`.

# Selection type

The selection type of [`UnarySelection`](@ref) can be any
`T <: AbstractContainer`. When queried for using the `selection_type` function,
will return the selection type of the given `sele`.

# Examples
```jldoctest
julia> sele = !rn"ALA"
UnarySelection ❯ ! "not" (Residue)
 └── FieldSelection › Residue.name = ALA
```
"""
mutable struct UnarySelection{M} <: AbstractSelection
    op::Function
    sele::AbstractSelection

    UnarySelection{M}(op::Function, sele::AbstractSelection) where {M <: AbstractStateMode} = begin
        new{M}(op, sele)
    end
end

state_mode_type(::UnarySelection{M}) where {M} = M
selection_type(sele::UnarySelection{M}) where {M} = selection_type(sele.sele)

# --- Unary Operations ---------------------------------------------------------
Base.:!(sele::AbstractSelection) = begin
    UnarySelection{state_mode_type(sele)}(!, sele)
end

# --- Select -------------------------------------------------------------------
function select(sele::UnarySelection{Stateless}, container::AbstractContainer)
    mask = select(sele.sele, container)
    mask = (sele.op)(mask)
end

function select(sele::UnarySelection{Stateful}, container::AbstractContainer)

    _selector = select(sele.sele, container)

    return function (state::State)
        mask = _selector(state)
        return (sele.op)(mask)
    end
end

# --- Show ---------------------------------------------------------------------

Base.show(io::IO, us::UnarySelection) = begin
    ProtoSyn.show(io, us)
end

function show(io::IO, us::UnarySelection{M}, levels::Opt{BitArray} = nothing) where {M}
    lead = ProtoSyn.get_lead(levels)
    dict = Dict((!) => "not")
    if levels === nothing
        levels = BitArray([])
    end
    println(io, lead*"UnarySelection ❯ $(us.op) \"$(dict[us.op])\" ($(selection_type(us)))")
    ProtoSyn.show(io, us.sele, vcat(levels, false))
end