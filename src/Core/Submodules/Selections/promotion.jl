export PromoteSelection
# Note: PromoteSelection is a BRANCH selection.

"""
    PromoteSelection(sele::AbstractSelection, ::Type{T}, op::Function) where {T <: AbstractContainer}

A [`PromoteSelection`](@ref) takes an input selection `sele` and outputs the
same result in a different [`Mask`](@ref) type, (depending on the operation `op`
requested for _upwards_ promotions). 

# State mode

The selection type of [`PromoteSelection`](@ref) can be either `Stateless` or
`Stateful`: it will automatically be set to the `StateMode` of the provided
`sele` on the constructor.

# Selection type

The selection type of [`PromoteSelection`](@ref) can be any
`T <: AbstractContainer`. When queried for using the `selection_type` function,
will return the selection type of the given `sele`.

!!! ukw "Note:"
    This selection does not have a short syntax version. However, the `promote`
    function is used to return [`PromoteSelection`](@ref) instances with a more
    user friendly syntax.

# Examples
```jldoctest
julia> sele = PromoteSelection(rn"ALA", Segment, all)
PromoteSelection ❯ From Residue to Segment
 └── FieldSelection › Residue.name = ALA
```
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

# --- Show ---------------------------------------------------------------------
function Base.show(io::IO, ps::PromoteSelection{M}, level_code::Opt{LevelCode} = nothing) where {M}
    lead = ProtoSyn.get_lead(level_code)
    if level_code === nothing
        level_code = LevelCode()
    end
    println(io, lead*"PromoteSelection ❯ From $(selection_type(ps.sele)) to $(ps.T)")
    Base.show(io, ps.sele, vcat(level_code, 4))
end

# ------------------------------------------------------------------------------

export promote

"""
    promote(sele::AbstractSelection, ::Type{T2}, [aggregator::Function = any]) where {T2 <: AbstractContainer}

Return a [`PromoteSelection`](@ref) instance for selection `sele`, promoting to
the requested type `T2 <: AbstractContainer`. If this is an _upwards_ promotion,
use the given `aggregator` function (default: `any`).

# Examples
```jldoctest
julia> ProtoSyn.promote(rn"ALA", Atom)
PromoteSelection ❯ From Residue to Atom
 └── FieldSelection › Residue.name = ALA
```
"""
function promote(sele::AbstractSelection, ::Type{T2}, aggregator::Function = any) where {T2 <: AbstractContainer}
    PromoteSelection(sele, T2, aggregator)
end


# --- Select -------------------------------------------------------------------
function select(sele::PromoteSelection{Stateless}, container::AbstractContainer)
    mask = (sele.sele)(container)
    return promote(mask, sele.T, container, sele.op)
end

function select(sele::PromoteSelection{Stateful}, container::AbstractContainer)

    _selector = select(sele.sele, container)

    return function (state::State)
        mask = _selector(state)
        return promote(mask, sele.T, container, sele.op)
    end
end


# --- Promotion tools for Masks ------------------------------------------------
"""
    promote(m1::Mask{T1}, m2::Mask{T2}, container::AbstractContainer) where {T1, T2}

Promote one of the 2 given `Masks` (`m1` and `m2`) to the lowest ranking common
type (Ex. `promote(m1::Mask{Residue}, m2::Mask{Atom}) => (Mask{Atom},
Mask{Atom)`).

# Examples
```jldoctest
julia> m1 = an"CB"(pose);

julia> m2 = rn"LEU"(pose);

julia> ProtoSyn.promote(m1, m2, pose.graph);
```
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

Cast a [`Mask`](@ref) of type `T1` to be of type `T2`, in the context of the
given `container`. If casting to a higher ranking type (_upwards_ promotion -
Ex. Atom -> Residue), a function `f` establishes the grouping operation (`any`
occurrence (by default) or `all` occurrences of lower ranking type must be
`true` to set the higher ranking entry to `true`.)

# Examples
```jldoctest
julia> m1 = an"CB"(pose);

julia> ProtoSyn.promote(m1, Residue, pose.graph);
```
"""
function promote(mask::Mask{T1}, ::Type{T2}, container::AbstractContainer, f::Function = any)::Mask where {T1 <: AbstractContainer, T2 <: AbstractContainer}

    new_mask = Mask{T2}(counter(T2)(container))

    if T1 > T2
        count = counter(T2)
        @assert T1 <= typeof(container) "Can't apply a promotion of $T1 to a single $(type(container)) instance."
        iter = T1 !== typeof(container) ? iterator(T1)(container) : [container]
        i = 1
        for (entry, item) in zip(mask, iter)
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