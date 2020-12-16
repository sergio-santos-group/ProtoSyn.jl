export PromoteSelection
# Note: PromoteSelection is a BRANCH selection.

"""
    PromoteSelection(sele::AbstractSelection, ::Type{T}, op::Function) where {T <: AbstractContainer}

A `PromoteSelection` takes an input selection `sele` and outputs the same result
in a different mask type, depending on the operation `op` requested. This is, in
essence, the same as calculating the `Mask` of a given `AbstractSelection` and
then using the function `promote` to cast the result to the desired `Mask` type
(Ex: `ProtoSyn.promote(mask, Type, container, f = op)`).  
The state mode of `PromoteSelection` `M` is set to the state mode of `sele`.

# Examples
```jldoctest
julia> sele = PromoteSelection(rn"ALA", Segment, all)
PromoteSelection{ProtoSyn.Stateless}(FieldSelection{ProtoSyn.Stateless,Residue}("ALA", :name, isequal), Segment, all)
```
```

## Note

Operation `op`is only employed in upwards promotions (Ex: `Atom` to `Residue`).
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


# --- Promotion tools ----------------------------------------------------------
"""
    promote(m1::Mask{T1}, m2::Mask{T2}, container::AbstractContainer) where {T1, T2}

Promote one of the 2 given `Masks` (`m1` and `m2`) to the lowest ranking common
type (Ex. `promote(m1::Mask{Residue}, m2::Mask{Atom}) => Mask{Atom}`).
"""
function promote(m1::Mask{T1}, m2::Mask{T2}, container::AbstractContainer) where {T1, T2}

    type1 = selection_type(m1)
    type2 = selection_type(m2)
    if type1 > type2
        m1 = promote(m1, type2, container)
    elseif type2 > type1
        m2 = promote(m2, type1, container)
    end

    m1, m2
end


"""
    promote(mask::Mask{T1}, ::Type{T2}, container::AbstractContainer, f::Function = any)::Mask where {T1 <: AbstractContainer, T2 <: AbstractContainer}

Cast a mask of type T1 to be of type T2, in the context of the given `container`.
If casting to a higher ranking type (Ex. Atom -> Residue), a function `f`
establishes the grouping operation (`any` occurrence (by default) or `all`
occurrences of lower ranking type must be `true` to set the higher ranking entry
to `true`.)

# Examples
```jldoctest
julia> promote(m1, Residue, container)
ProtoSyn.Mask{Residue}(Bool[1, 1, 1, 1, 1])
```
"""
function promote(mask::Mask{T1}, ::Type{T2}, container::AbstractContainer, f::Function = any)::Mask where {T1 <: AbstractContainer, T2 <: AbstractContainer}

    new_mask = Mask{T2}(counter(T2)(container))

    if T1 > T2
        count = counter(T2)
        iter = iterator(T1)
        i = 1
        for (entry, item) in zip(mask, iter(container))
            n = count(item)
            new_mask[i:(i+n-1)] = entry
            i += n
        end
    else
        count = counter(T1)
        iter = iterator(T2)
        i = 1
        for (index, item) in enumerate(iter(container))
            n = count(item)
            new_mask[index] = f(view(mask.content, i:(i+n-1)))
            i += n
        end
    end

    return new_mask
end