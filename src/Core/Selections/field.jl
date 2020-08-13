export FieldSelection
# Note: FieldSelection is a LEAF selection.

"""
    FieldSelection{M, T} <: AbstractSelection

A `FieldSelection` takes an input `pattern` (as a `String`) and a `field` (as a
`Symbol`) and outputs a `Mask` (of type `T <: AbstractContainer`) containing all
instances of said type in the given `container` whose `field` matches the `pattern`.

    FieldSelection{T}(pattern::String, field::Symbol) where {T <: AbstractContainer}
    
The given `pattern` is to be considered as a Regular Expression (`Regex`).
The state mode of DistanceSelection `M` is forced to be Stateless.
The short syntax macros are: an"..." = Atom name; as"..." = Atom symbol; rn"..."
= Residue name; sn"..." = Segment name

# Examples
```jldoctest
julia> sele = FieldSelection{Segment}("UNK", :name)
FieldSelection{ProtoSyn.Stateless,Segment}(true, r"UNK", :name)

julia> sele = FieldSelection{Residue}("ALA", :name)
FieldSelection{ProtoSyn.Stateless,Residue}(true, r"ALA", :name)

julia> sele = FieldSelection{Residue}("GL*", :name)
FieldSelection{ProtoSyn.Stateless,Residue}(true, r"GL*", :name)

julia> sele = FieldSelection{Atom}("CA", :name)
FieldSelection{ProtoSyn.Stateless,Atom}(true, r"CA", :name)

julia> sele = FieldSelection{Atom}("C", :symbol)
FieldSelection{ProtoSyn.Stateless,Atom}(true, r"C", :symbol)
```

# Short Syntax
```jldoctest
julia> sele = sn"UNK"
FieldSelection{ProtoSyn.Stateless,Segment}(true, r"UNK", :name)

julia> sele = rn"ALA"
FieldSelection{ProtoSyn.Stateless,Residue}(true, r"ALA", :name)

julia> sele = rn"GL*"
FieldSelection{ProtoSyn.Stateless,Residue}(true, r"GL*", :name)

julia> sele = an"CA"
FieldSelection{ProtoSyn.Stateless,Atom}(true, r"CA", :name)

julia> sele = as"C"
FieldSelection{ProtoSyn.Stateless,Atom}(true, r"C", :symbol)
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