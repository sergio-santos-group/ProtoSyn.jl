export SerialSelection
# Note: SerialSelection is a LEAF selection.

"""
    SerialSelection{T}(serial::Int, field::Symbol) where {T <: AbstractContainer}
    

A [`SerialSelection`](@ref) selects instances based on `:id` and `:index`. It
takes an input `serial` (as an `Int`) and a `field` (as a `Symbol`) and outputs
a [`Mask`](@ref) (of type `T <: AbstractContainer`) containing all instances of
said type in the given `container` whose `field` matches the `serial` number
given marked as `true`. This selection works similarly to
[`FieldSelection`](@ref), but is especialized in dealing with number variables.

# State mode

The state mode of [`SerialSelection`](@ref) `M` is forced to be `Stateless`.

# Selection type

The selection type of [`SerialSelection`](@ref) can be any `T <: AbstractContainer`.

# Short syntax

* sid"..." = Segment ID;
* rid"..." = Residue ID;
* aid"..." = Atom ID;
* six"..." = Segment Index;
* rix"..." = Residue Index;
* aix"..." = Atom Index;

# Examples
```jldoctest
julia> sele = SerialSelection{Segment}(1, :id)
SerialSelection › Segment.id = 1

julia> sele = SerialSelection{Atom}(2, :index)
SerialSelection › Atom.index = 2

julia> sele = sid"1"
SerialSelection › Segment.id = 1

julia> sele = aix"2"
SerialSelection › Atom.index = 2
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

    iter = typeof(container) > T ? iterator(T)(container) : [container]
    for item in iter
        if sele.serial == getproperty(item, sele.field)
            mask[item.index] = true
        end
    end

    return mask
end

# --- Show ---------------------------------------------------------------------
function Base.show(io::IO, ss::SerialSelection{M, T}, level_code::Opt{LevelCode} = nothing) where {M, T}
    lead = ProtoSyn.get_lead(level_code)
    if level_code === nothing
        level_code = LevelCode()
    end
    println(io, lead*"SerialSelection › $(T).$(ss.field) = $(ss.serial)")
end

# --- Help ---------------------------------------------------------------------
SerialSelection(serial::Int, field::Symbol) = begin
    @error """
    SerialSelection requires specification of selection type: Atom, Residue or Segment. Some examples:
     • SerialSelection{Atom}($serial, :$field)
     • SerialSelection{Residue}($serial, :$field)
     • SerialSelection{Segment}($serial, :$field)
    """
end

# ------------------------------------------------------------------------------

export RangeSelection
# Note: RangeSelection is a LEAF selection.

"""
    RangeSelection{T}(range::UnitRange{Int}, field::Symbol) where {T <: AbstractContainer}

A [`RangeSelection`](@ref) takes an input `range` (as an `UnitRange{Int}`) and a `field`
(as a `Symbol`) and outputs a [`Mask`](@ref) (of type `T <: AbstractContainer`)
containing all instances of said type in the given `container` whose `field`
matches is in the `range` given. This selection works simillarly to
[`FieldSelection`](@ref), but is especialized in dealing with numbers.

!!! ukw "Note:"
    The [`RangeSelection`](@ref) is inclusive, meaning the `:id` of `:index`
    given in the selection will also be included in the selected [`Mask`](@ref).

# State mode

The state mode of [`RangeSelection`](@ref) `M` is forced to be `Stateless`.

# Selection type

The selection type of [`RangeSelection`](@ref) can be any `T <: AbstractContainer`.

# Short syntax
* sid"..." = Segment ID
* rid"..." = Residue ID
* aid"..." = Atom ID
* six"..." = Segment Index;
* rix"..." = Residue Index;
* aix"..." = Atom Index;

# Examples
```jldoctest
julia> sele = RangeSelection{Segment}(1:4, :id)
RangeSelection › Segment.id between 1 and 4

julia> sele = RangeSelection{Atom}(2:10, :index)
RangeSelection › Atom.index between 2 and 10

julia> sele = sid"1:4"
RangeSelection › Segment.id between 1 and 4

julia> sele = aix"2:10"
RangeSelection › Atom.index between 2 and 10
```

!!! ukw "Note:"
    This selection assumes that all `Abstractcontainer` instances are ordered
    (i.e: `aid"1:10"` will select atoms 1, 2, 3, 4, 5, 6, 7, 8, 9 and 10).
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

    s1 = s[1]
    s2 = s[2]
    parse(Int, s2) < parse(Int, s1) && begin
        s1 = s[2]
        s2 = s[1]
        @warn "Range should be ordered (\"x:y\" where x < y). The provided range (\"$string\") was inverted (to \"$s1:$s2\")."
    end
    return parse(Int, s1):parse(Int, s2)
end

export @aid_str, @aix_str, @sid_str, @six_str, @rid_str, @rix_str
macro sid_str(p); occursin(":", p) && return RangeSelection{Segment}(parse(UnitRange{Int}, p), :id); SerialSelection{Segment}(parse(Int, p), :id); end
macro six_str(p); occursin(":", p) && return RangeSelection{Segment}(parse(UnitRange{Int}, p), :index); SerialSelection{Segment}(parse(Int, p), :index); end
macro rid_str(p); occursin(":", p) && return RangeSelection{Residue}(parse(UnitRange{Int}, p), :id); SerialSelection{Residue}(parse(Int, p), :id); end
macro rix_str(p); occursin(":", p) && return RangeSelection{Residue}(parse(UnitRange{Int}, p), :index); SerialSelection{Residue}(parse(Int, p), :index); end
macro aid_str(p); occursin(":", p) && return RangeSelection{Atom}(parse(UnitRange{Int}, p), :id); SerialSelection{Atom}(parse(Int, p), :id); end
macro aix_str(p); occursin(":", p) && return RangeSelection{Atom}(parse(UnitRange{Int}, p), :index); SerialSelection{Atom}(parse(Int, p), :index); end

# --- Show ---------------------------------------------------------------------
function Base.show(io::IO, rs::RangeSelection{M, T}, level_code::Opt{LevelCode} = nothing) where {M, T}
    lead = ProtoSyn.get_lead(level_code)
    if level_code === nothing
        level_code = LevelCode()
    end
    if typeof(rs.range) === Int
        println(io, lead*"RangeSelection › $(T).$(rs.field) larger than $(rs.range)")
    else
        println(io, lead*"RangeSelection › $(T).$(rs.field) between $(rs.range.start) and $(rs.range.stop)")
    end
end