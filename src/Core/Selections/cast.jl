export PromoteSelection
# Note: PromoteSelection is a BRANCH selection.

"""
    PromoteSelection{M} <: AbstractSelection

A `PromoteSelection` takes an input selection `sele` and outputs the same result in a different
mask type, depending on the operation `op` requested. This is, in essence, the same
as calculating the `Mask` of a given `AbstractSelection` and then using the function `promote`
to cast the result to the desired `Mask` type (Ex: `promote(mask, Type, container, f = all)`)

    PromoteSelection(op::Function, sele::AbstractSelection, ::Type{T}) where {M <: AbstractStateMode, T <: AbstractContainer}
    
The state mode of `PromoteSelection` `M` is set to the state mode of `sele`.

# Examples
```jldoctest
julia> sele = PromoteSelection(all, rn"ALA", Segment)
PromoteSelection{ProtoSyn.Stateless}(true, any, FieldSelection{ProtoSyn.Stateless,Residue}(true, r"ALA", :name), Segment)

julia> sele = PromoteSelection(any, an"CA", Residue)
PromoteSelection{ProtoSyn.Stateless}(true, any, FieldSelection{ProtoSyn.Stateless,Atom}(false, r"CA", :name), Residue)

julia> sele = PromoteSelection(all, rn"ALA", Atom)
PromoteSelection{ProtoSyn.Stateless}(true, all, FieldSelection{ProtoSyn.Stateless,Residue}(true, r"ALA", :name), Atom)
```

# Short Syntax
```jldoctest
julia> all(rn"ALA", Segment)
PromoteSelection{ProtoSyn.Stateless}(true, any, FieldSelection{ProtoSyn.Stateless,Residue}(true, r"ALA", :name), Segment)

julia> any(an"CA", Residue)
PromoteSelection{ProtoSyn.Stateless}(true, any, FieldSelection{ProtoSyn.Stateless,Atom}(false, r"CA", :name), Residue)

julia> cast(rn"ALA", Atom)
PromoteSelection{ProtoSyn.Stateless}(true, all, FieldSelection{ProtoSyn.Stateless,Residue}(true, r"ALA", :name), Atom)
```

## Note

When using `any` or `all`, only upwards promotions are possible (Ex: `Atom` to
`Residue`). Similarly, when using `cast` only downwards promotions are possible
(Ex: `Residue` to `Atom`).
"""
mutable struct PromoteSelection{M} <: AbstractSelection
    sele::AbstractSelection
    T::Type{<: AbstractContainer}
    op::Function

    PromoteSelection(sele::AbstractSelection, ::Type{T}, op::Function) where {T <: AbstractContainer} = begin
        new{state_mode_type(sele)}(sele, T, op)
    end
end

state_mode_type(::PromoteSelection{M}) where {M} = M
selection_type(sele::PromoteSelection{M})  where {M} = sele.T


export promote
function promote(sele::AbstractSelection, ::Type{T2}, aggregator::Function = any) where {T2 <: AbstractContainer}
    PromoteSelection(sele, T2, aggregator)
end


# --- Select -------------------------------------------------------------------
function select(sele::PromoteSelection{Stateless}, container::AbstractContainer)
    mask = select(sele.sele, container)
    return promote(mask, sele.T, container, sele.op)
end

function select(sele::PromoteSelection{Stateful}, container::AbstractContainer)

    _selector = select(sele.sele, container)

    return function (state::State)
        mask = _selector(state)
        return promote(mask, sele.T, container, sele.op)
    end
end