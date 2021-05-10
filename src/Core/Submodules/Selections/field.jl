export FieldSelection
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
mutable struct FieldSelection{M, T} <: AbstractSelection
    pattern::Union{Regex, String}
    field::Symbol
    op::Function
    
    FieldSelection{T}(pattern::String, field::Symbol; is_regex = false) where {T <: AbstractContainer} = begin
        if is_regex
            new{Stateless, T}(Regex(pattern), field, occursin)
        else
            new{Stateless, T}(pattern, field, isequal)
        end
    end
end

state_mode_type(::FieldSelection{M, T}) where {M, T} = M
selection_type(::FieldSelection{M, T}) where {M, T} = T

# --- Select -------------------------------------------------------------------
function select(sele::FieldSelection{Stateless, T}, container::AbstractContainer) where {T <: AbstractContainer}

    n_items = counter(T)(container)
    mask = Mask{T}(n_items)

    for (index, item) in enumerate(iterator(T)(container))
        if sele.op(sele.pattern, getproperty(item, sele.field))
            mask[index] = true
        end
    end
    return mask
end

# --- Short Syntax -------------------------------------------------------------
export @sn_str, @rn_str, @an_str, @as_str
parse_flags(flags) = isempty(flags) ? false : 'r' in flags[1]
macro sn_str(p, flags...); FieldSelection{Segment}(p, :name; is_regex = parse_flags(flags)); end
macro rn_str(p, flags...); FieldSelection{Residue}(p, :name; is_regex = parse_flags(flags)); end
macro an_str(p, flags...); FieldSelection{Atom}(p, :name; is_regex = parse_flags(flags)); end
macro as_str(p, flags...); FieldSelection{Atom}(p, :symbol; is_regex = parse_flags(flags)); end

# --- Show ---------------------------------------------------------------------

Base.show(io::IO, fs::FieldSelection) = begin
    ProtoSyn.show(io, fs)
end

function show(io::IO, fs::FieldSelection{M, T}, levels::Opt{BitArray} = nothing) where {M, T}
    lead = ProtoSyn.get_lead(levels)
    if levels === nothing
        levels = BitArray([])
    end
    println(io, lead*"FieldSelection › $(T).$(fs.field) = $(fs.pattern)")
end