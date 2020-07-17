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
    
The state mode of DistanceSelection `M` is forced to be Stateless.
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
    is_exit_node::Bool
    serial::Int
    field::Symbol
    
    SerialSelection{T}(serial::Int, field::Symbol) where {T <: AbstractContainer} = begin
        new{Stateless, T}(true, serial, field)
    end
end

state_mode_type(::SerialSelection{M, T}) where {M, T} = M
selection_type(::SerialSelection{M, T}) where {M, T} = T

# --- Select -------------------------------------------------------------------
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

# --- Short Syntax -------------------------------------------------------------
export @aid_str, @aix_str, @sid_str, @six_str, @rid_str, @rix_str
macro sid_str(p); SerialSelection{Segment}(parse(Int, p), :id); end
macro six_str(p); SerialSelection{Segment}(parse(Int, p), :index); end
macro rid_str(p); SerialSelection{Residue}(parse(Int, p), :id); end
macro rix_str(p); SerialSelection{Residue}(parse(Int, p), :index); end
macro aid_str(p); SerialSelection{Atom}(parse(Int, p), :id); end
macro aix_str(p); SerialSelection{Atom}(parse(Int, p), :index); end

# ---

export MaxSerialSelection
mutable struct MaxSerialSelection{M, T} <: AbstractSelection
    is_exit_node::Bool
    field::Symbol
    
    MaxSerialSelection{T}(field::Symbol) where {T <: AbstractContainer} = begin
        new{Stateless, T}(true, field)
    end
end

state_mode_type(::MaxSerialSelection{M, T}) where {M, T} = M
selection_type(::MaxSerialSelection{M, T}) where {M, T} = T

# --- Select -------------------------------------------------------------------
function select(sele::MaxSerialSelection{Stateless, T}, container::AbstractContainer) where {T <: AbstractContainer}

    n_items = counter(T)(container)
    mask = Mask{T}(n_items)

    max_serial = 0
    for item in iterator(T)(container)
        prop = getproperty(item, sele.field)
        println(prop)
        if prop > max_serial
            mask = Mask{T}(n_items)
            mask[item.index] = true
            max_serial = prop
        elseif prop == max_serial
            mask[item.index] = true
        end
    end
    return mask
end
