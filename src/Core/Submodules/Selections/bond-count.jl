export BondCountSelection
# Note: FieldSelection is a LEAF selection.

"""
    FieldSelection{T}(pattern::String, field::Symbol, [is_regex::Bool = false]) where {T <: AbstractContainer}

A [`FieldSelection`](@ref) takes an input `pattern` (as a `String`) and a `field` (as a
`Symbol`) and outputs a `Mask` (of type `T <: AbstractContainer`) containing all
instances of said type in the given `container` whose `field` matches the `pattern`.

!!! ukw "Note:"
    The given `pattern` can be considered as a Regular Expression (`Regex`), if
    `is_regex` flag is set to `true`. Optinally, when using a short syntax,
    appending an "r" flag at the end of the expression also sets `is_regex` to
    `true`. Check the examples bellow.

# State mode

The state mode of [`FieldSelection`](@ref) `M` is forced to be `Stateless`.

# Selection type

The selection type of [`FieldSelection`](@ref) can be any `T <: AbstractContainer`.

# Short syntax
* an"..." = Atom name
* as"..." = Atom symbol
* rn"..." = Residue name
* sn"..." = Segment name

# Examples
```jldoctest
julia> sele = FieldSelection{Atom}("C", :symbol)
FieldSelection › Atom.symbol = C

julia> sele = FieldSelection{Residue}("AL*", :name, is_regex = true)
FieldSelection › Residue.name = r"AL*"

julia> sele = as"C"
FieldSelection › Atom.symbol = C

julia> sele = rn"AL*"r
FieldSelection › Residue.name = r"AL*"
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