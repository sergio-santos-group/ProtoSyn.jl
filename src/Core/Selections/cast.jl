export CastSelection
# Note: CastSelection is a BRANCH selection.

"""
    CastSelection{M} <: AbstractSelection

A `CastSelection` takes an input selection `sele` and outputs the same result in a different
mask type, depending on the operation `op` requested. This is, in essence, the same
as calculating the `Mask` of a given `AbstractSelection` and then using the function `promote`
to cast the result to the desired `Mask` type (Ex: `promote(mask, Type, container, f = all)`)

    CastSelection(op::Function, sele::AbstractSelection, ::Type{T}) where {M <: AbstractStateMode, T <: AbstractContainer}
    
The state mode of `CastSelection` `M` is set to the state mode of `sele`.

# Examples
```jldoctest
julia> sele = CastSelection(all, rn"ALA", Segment)
CastSelection{ProtoSyn.Stateless}(true, any, FieldSelection{ProtoSyn.Stateless,Residue}(true, r"ALA", :name), Segment)

julia> sele = CastSelection(any, an"CA", Residue)
CastSelection{ProtoSyn.Stateless}(true, any, FieldSelection{ProtoSyn.Stateless,Atom}(false, r"CA", :name), Residue)

julia> sele = CastSelection(all, rn"ALA", Atom)
CastSelection{ProtoSyn.Stateless}(true, all, FieldSelection{ProtoSyn.Stateless,Residue}(true, r"ALA", :name), Atom)
```

# Short Syntax
```jldoctest
julia> all(rn"ALA", Segment)
CastSelection{ProtoSyn.Stateless}(true, any, FieldSelection{ProtoSyn.Stateless,Residue}(true, r"ALA", :name), Segment)

julia> any(an"CA", Residue)
CastSelection{ProtoSyn.Stateless}(true, any, FieldSelection{ProtoSyn.Stateless,Atom}(false, r"CA", :name), Residue)

julia> cast(rn"ALA", Atom)
CastSelection{ProtoSyn.Stateless}(true, all, FieldSelection{ProtoSyn.Stateless,Residue}(true, r"ALA", :name), Atom)
```

## Note

When using `any` or `all`, only upwards promotions are possible (Ex: `Atom` to
`Residue`). Similarly, when using `cast` only downwards promotions are possible
(Ex: `Residue` to `Atom`).
"""
mutable struct CastSelection{M} <: AbstractSelection
    is_exit_node::Bool
    op::Function
    sele::AbstractSelection
    T::Type{<: AbstractContainer}

    CastSelection(op::Function, sele::AbstractSelection, ::Type{T}) where {T <: AbstractContainer} = begin
        new{state_mode_type(sele)}(true, op, sele, T)
    end
end

state_mode_type(::CastSelection{M}) where {M} = M
selection_type(sele::CastSelection{M})  where {M} = sele.T

# --- Short Syntax -------------------------------------------------------------
Base.any(sele::AbstractSelection, ::Type{T2}) where {T2 <: AbstractContainer} = begin
    T1 = selection_type(sele)
    @assert T1 < T2 "'any' function casts to higher order AbstractContainers only"
    sele.is_exit_node = false
    CastSelection(any, sele, T2)
end


Base.all(sele::AbstractSelection, ::Type{T2}) where {T2 <: AbstractContainer} = begin
    T1 = selection_type(sele)
    @assert T1 < T2 "'any' function casts to higher order AbstractContainers only"
    sele.is_exit_node = false
    CastSelection(all, sele, T2)
end


export cast
function cast(sele::AbstractSelection, ::Type{T2}) where {T2 <: AbstractContainer}
    T1 = selection_type(sele)
    @assert T1 > T2 "'cast' function casts to lower order AbstractContainers only"
    sele.is_exit_node = false
    CastSelection(all, sele, T2)
end


# --- Select -------------------------------------------------------------------
function select(sele::CastSelection{Stateless}, container::AbstractContainer)
    mask = select(sele.sele, container)
    return promote(mask, sele.T, container, sele.op)
end

function select(sele::CastSelection{Stateful}, container::AbstractContainer)

    _selector = select(sele.sele, container)

    return function (state::State)
        mask = _selector(state)
        return promote(mask, sele.T, container, sele.op)
    end
end