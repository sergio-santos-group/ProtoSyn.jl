export UnarySelection
# Note: UnarySelection is a BRANCH selection.

"""
    UnarySelection{M} <: AbstractSelection

A `UnarySelection` applies an operation `op` to the given `AbstractSelection`
`sele`, such as `!`.

# Examples
```jldoctest
julia> sele = !rn"ALA"
UnarySelection{ProtoSyn.Stateless}(true, !, FieldSelection{ProtoSyn.Stateless,Residue}(false, r"ALA", :name))
```
"""
mutable struct UnarySelection{M} <: AbstractSelection
    is_exit_node::Bool
    op::Function
    sele::AbstractSelection

    UnarySelection{M}(op::Function, sele::AbstractSelection) where {M <: AbstractStateMode} = begin
        new{M}(true, op, sele)
    end
end

state_mode_type(::UnarySelection{M}) where {M} = M
selection_type(sele::UnarySelection{M}) where {M} = selection_type(sele.sele)

# --- Unary Operations ---------------------------------------------------------
Base.:!(sele::AbstractSelection) = begin
    sele.is_exit_node = false
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