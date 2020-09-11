export SerialSelection
# Note: SerialSelection is a LEAF selection.

"""
    SerialSelection{M, T} <: AbstractSelection

A `SerialSelection` takes an input `serial` (as an `Int`) and a `field` (as a
`Symbol`) and outputs a `Mask` (of type `T <: AbstractContainer`) containing all
instances of said type in the given `container` whose `field` matches the
`serial` number given. This selection works simillarly to FieldSelection, except
it is especialized in dealing with numbers.

    SerialSelection{T}(serial::Int, field::Symbol) where {T <: AbstractContainer}
    
The state mode of SerialSelection `M` is forced to be Stateless.
The short syntax macros are: sid"..." = Segment ID; six"..." = Segment Index;
rid"..." = Residue ID; rix"..." = Residue Index;
aid"..." = Atom ID; aix"..." = Atom Index;

# Examples
```jldoctest
julia> sele = SerialSelection{Segment}(1, :id)
SerialSelection{ProtoSyn.Stateless,Segment}(true, 1, :id)

julia> sele = SerialSelection{Atom}(2, :index)
SerialSelection{ProtoSyn.Stateless,Atom}(true, 2, :index)
```

# Short Syntax
```jldoctest
julia> sele = sid"1"
SerialSelection{ProtoSyn.Stateless,Segment}(true, 1, :id)

julia> sele = aix"2"
SerialSelection{ProtoSyn.Stateless,Atom}(true, 2, :index)
```
"""
mutable struct SerialSelection{M, T} <: AbstractSelection
    serial::Int
    field::Symbol
    
    SerialSelection{T}(serial::Int, field::Symbol) where {T <: AbstractContainer} = begin
        new{Stateless, T}(serial, field)
    end
end

state_mode_type(::SerialSelection{M, T}) where {M, T} = M
selection_type(::SerialSelection{M, T}) where {M, T} = T

# --- SerialSelection Select ---------------------------------------------------
function select(sele::SerialSelection{Stateless, T}, container::AbstractContainer) where {T <: AbstractContainer}

    n_items = counter(T)(container)
    mask = Mask{T}(n_items)

    for item in iterator(T)(container)
        if sele.serial == getproperty(item, sele.field)
            mask[item.index] = true
        end
    end
    return mask
end

# ------------------------------------------------------------------------------

export RangeSelection
# Note: RangeSelection is a LEAF selection.

"""
    RangeSelection{M, T} <: AbstractSelection

A `RangeSelection` takes an input `range` (as an `UnitRange{Int}`) and a `field`
(as a `Symbol`) and outputs a `Mask` (of type `T <: AbstractContainer`)
containing all instances of said type in the given `container` whose `field`
matches is in the `range` given. This selection works simillarly to
FieldSelection, except it is especialized in dealing with numbers.

    RangeSelection{T}(range::UnitRange{Int}, field::Symbol) where {T <: AbstractContainer}
    
The state mode of RangeSelection `M` is forced to be Stateless.
# The short syntax macros are: sid"..." = Segment ID; six"..." = Segment Index;
# rid"..." = Residue ID; rix"..." = Residue Index;
# aid"..." = Atom ID; aix"..." = Atom Index;

# Examples
```jldoctest
julia> sele = SerialSelection{Segment}(1, :id)
SerialSelection{ProtoSyn.Stateless,Segment}(true, 1, :id)

julia> sele = SerialSelection{Atom}(2, :index)
SerialSelection{ProtoSyn.Stateless,Atom}(true, 2, :index)
```

# Short Syntax
```jldoctest
julia> sele = sid"1"
SerialSelection{ProtoSyn.Stateless,Segment}(true, 1, :id)

julia> sele = aix"2"
SerialSelection{ProtoSyn.Stateless,Atom}(true, 2, :index)
```
"""
mutable struct RangeSelection{M, T} <: AbstractSelection
    range::Union{UnitRange{Int}, Int}
    field::Symbol
    
    RangeSelection{T}(range::Union{UnitRange{Int}, Int}, field::Symbol) where {T <: AbstractContainer} = begin
        new{Stateless, T}(range, field)
    end
end

state_mode_type(::RangeSelection{M, T}) where {M, T} = M
selection_type(::RangeSelection{M, T}) where {M, T} = T

# --- RangeSelection Select ----------------------------------------------------
function select(sele::RangeSelection{Stateless, T}, container::AbstractContainer) where {T <: AbstractContainer}

    n_items = counter(T)(container)
    mask = Mask{T}(n_items)

    if typeof(sele.range) == Int
        range = sele.range:n_items
    else
        range = sele.range
    end

    for item in iterator(T)(container)
        if getproperty(item, sele.field) in range
            mask[item.index] = true
        end
    end
    return mask
end

# --- Short Syntax -------------------------------------------------------------

Base.parse(::Type{UnitRange{Int}}, string::String; base = 2) = begin
    s = split(string, ":")
    if s[2] == "end"
        return parse(Int, s[1])
    end
    return parse(Int, s[1]):parse(Int, s[2])
end

export @aid_str, @aix_str, @sid_str, @six_str, @rid_str, @rix_str
macro sid_str(p); occursin(":", p) && return RangeSelection{Segment}(parse(UnitRange{Int}, p), :id); SerialSelection{Segment}(parse(Int, p), :id); end
macro six_str(p); occursin(":", p) && return RangeSelection{Segment}(parse(UnitRange{Int}, p), :index); SerialSelection{Segment}(parse(Int, p), :index); end
macro rid_str(p); occursin(":", p) && return RangeSelection{Residue}(parse(UnitRange{Int}, p), :id); SerialSelection{Residue}(parse(Int, p), :id); end
macro rix_str(p); occursin(":", p) && return RangeSelection{Residue}(parse(UnitRange{Int}, p), :index); SerialSelection{Residue}(parse(Int, p), :index); end
macro aid_str(p); occursin(":", p) && return RangeSelection{Atom}(parse(UnitRange{Int}, p), :id); SerialSelection{Atom}(parse(Int, p), :id); end
macro aix_str(p); occursin(":", p) && return RangeSelection{Atom}(parse(UnitRange{Int}, p), :index); SerialSelection{Atom}(parse(Int, p), :index); end