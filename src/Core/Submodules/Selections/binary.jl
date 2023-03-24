export BinarySelection
# Note: BinarySelection is a BRANCH selection.

"""
    BinarySelection(op::Function, left::L, right::R) where {L <: AbstractSelection, R <: AbstractSelection}

A [`BinarySelection`](@ref) merges two selections using different operators
`op`, such as `and` and `or`. Return a new [`BinarySelection`](@ref) that
combines both `left` and `right` AbstractSelections, using the defined operator
`op`.

# State mode

If the defined selections have the same state mode, the resulting mask will be
of that state mode. If the defined selections have different state modes (Ex:
`Stateless` and `Stateful`) the resulting selection will have a `Stateful` state
mode.

# Selection type

If the defined have the same selection type, the resulting mask will be of that
type. If the defined selections have different selection types (Ex: `Atom`
and `Residue`), the resulting mask will be promoted to the lowest ranking type
(Ex: `Atom`).

# See also
[Promotion](@ref)

# Short syntax
* ... & ...
* ... | ...

!!! ukw "Note:"
    [`BinarySelection`](@ref) instances are _left-dominant_, meaning that a
    grouping of logical operators such as `rn"ALA" & rn"LEU" | an"CA"` will
    first resolve the `rn"ALA" & rn"LEU"` part (which should return an all-false
    [`Mask`](@ref)) and only then combine this [`Mask`](@ref) with the `| an"CA"`
    [`BinarySelection`](@ref), thus essentially selecting only the _CA_ atoms of
    the `AbstractContainer` its applied to.

    However, **selections respect to parenthesis grouping**, meaning `rn"ALA" &
    (rn"LEU" | an"CA")` will first resolve `rn"LEU" | an"CA"` (which should
    return a [`Mask`](@ref) selecting all atoms of all _LEU_ residues plus the
    _CA_ atoms of all other residues) and only then combine this [`Mask`](@ref)
    with the `rn"ALA" &` [`BinarySelection`](@ref), thus essentially selecting
    only the _CA_ atoms of any _ALA_ residues of the `AbstractContainer` its
    applied to.

# Examples
```jldoctest
julia> sele = BinarySelection(&, rn"ALA", an"CA")
BinarySelection ❯  & "and" (Atom)
 ├── FieldSelection › Residue.name = ALA
 └── FieldSelection › Atom.name = CA

julia> rn"ALA" & an"CA"
BinarySelection ❯  & "and" (Atom)
 ├── FieldSelection › Residue.name = ALA
 └── FieldSelection › Atom.name = CA
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

# Case matching with nothing
Base.:&(::Nothing, sele::AbstractSelection) = sele
Base.:&(sele::AbstractSelection, ::Nothing) = sele
Base.:|(::Nothing, sele::AbstractSelection) = sele
Base.:|(sele::AbstractSelection, ::Nothing) = sele

Base.:&(::Nothing, ::Nothing) = nothing
Base.:|(::Nothing, ::Nothing) = nothing

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

# --- Show ---------------------------------------------------------------------
function Base.show(io::IO, bs::BinarySelection{M, T}, level_code::Opt{LevelCode} = nothing) where {M, T}
    lead = ProtoSyn.get_lead(level_code)
    dict = Dict((|) => "or", (&) => "and")
    if level_code === nothing
        level_code = LevelCode()
    end
    if length(level_code.levels) > 5
        println(io, lead*"(...)")
    else
        println(io, lead*"BinarySelection ❯  $(bs.op) \"$(dict[bs.op])\" ($(selection_type(bs)))")
        Base.show(io, bs.left, vcat(level_code, 3))
        Base.show(io, bs.right, vcat(level_code, 4))
    end
end