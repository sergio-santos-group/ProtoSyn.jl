export BinarySelection
# Note: BinarySelection is a BRANCH selection.

"""
    BinarySelection{LM, RM} <: AbstractSelection

A `BinarySelection` merges two selections using different operators, such as `and`
and `or`.

    BinarySelection(op::Function, left::L, right::R) where {L <: AbstractSelection, R <: AbstractSelection}
    
Construct a new `BinarySelection` that combines both `left` and `right` AbstractSelections,
using the defined operator `op`. If the defined selections have different resulting
mask types (Ex: `Atom` and `Residue`), the resulting mask will be promoted to the lowest
ranking type (Ex: `Atom`). If the defined selections have different state modes
(Ex: `Stateless` and `Stateful`) the resulting selection will have a `Stateful` state mode.

# Examples
```jldoctest
julia> sele = BinarySelection(&, rn"ALA", an"CA")
BinarySelection{ProtoSyn.Stateless,ProtoSyn.Stateless}(true, &, FieldSelection{ProtoSyn.Stateless,Atom}(true, r"CA", :name), FieldSelection{ProtoSyn.Stateless,Residue}(true, r"ALA", :name))
```

# Short Syntax
```jldoctest
julia> rn"ALA" & an"CA"
BinarySelection{ProtoSyn.Stateless,ProtoSyn.Stateless}(true, &, FieldSelection{ProtoSyn.Stateless,Atom}(true, r"CA", :name), FieldSelection{ProtoSyn.Stateless,Residue}(true, r"ALA", :name))
```
"""
mutable struct BinarySelection{LM, RM} <: AbstractSelection
    op::Function
    left::AbstractSelection
    right::AbstractSelection

    BinarySelection(op::Function, left::L, right::R) where {L <: AbstractSelection, R <: AbstractSelection} = begin
        
        ML = state_mode_type(left)
        MR = state_mode_type(right)

        new{ML, MR}(op, left, right)
    end
end


# --- State Mode Rules ---------------------------------------------------------
state_mode_rule(::Type{Stateless}, ::Type{Stateless}) = Stateless
state_mode_rule(::Type{Stateful},  ::Type{Stateless}) = Stateful
state_mode_rule(::Type{Stateless}, ::Type{Stateful})  = Stateful
state_mode_rule(::Type{Stateful},  ::Type{Stateful})  = Stateful

state_mode_type(::BinarySelection{LM, RM}) where {LM, RM} = state_mode_rule(LM, RM)
selection_type(sele::BinarySelection) = min(selection_type(sele.left), selection_type(sele.right))

# --- Short Syntax -------------------------------------------------------------
Base.:&(l::AbstractSelection, r::AbstractSelection) = BinarySelection(&, l, r)
Base.:|(l::AbstractSelection, r::AbstractSelection) = BinarySelection(|, l, r)


# --- Select -------------------------------------------------------------------
select(sele::BinarySelection{Stateless, Stateless}, container::AbstractContainer) = begin
    
    left_mask  = select(sele.left,  container)
    right_mask = select(sele.right, container)

    l_mask, r_mask = promote(left_mask, right_mask, container)
    
    return sele.op(l_mask, r_mask)
end

select(sele::BinarySelection{Stateful, Stateless}, container::AbstractContainer) = select(sele, container, sele.right, sele.left)

select(sele::BinarySelection{Stateless, Stateful}, container::AbstractContainer) = select(sele, container, sele.left, sele.right)

select(sele::BinarySelection{T1, T2}, container::AbstractContainer, sless::AbstractSelection, sful::AbstractSelection) where {T1, T2} = begin
    
    sless_mask  = select(sless,  container)
    sful_selector = select(sful, container)

    return function (state::State)
        sful_mask = sful_selector(state)
        slessmask, sfulmask = promote(sless_mask, sful_mask, container)
        return sele.op(slessmask, sfulmask)
    end
end

select(sele::BinarySelection{Stateful, Stateful}, container::AbstractContainer) = begin
    
    left_selector  = select(sele.left,  container)
    right_selector = select(sele.right, container)

    return function (state::State)
        lmask = left_selector(state)
        rmask = right_selector(state)
        l_mask, r_mask = promote(lmask, rmask, container)
        return sele.op(l_mask, r_mask)
    end
end