export RandomSelection
# Note: RandomSelection is a LEAF selection.

"""
    RandomSelection{M, T} <: AbstractSelection

A `RandomSelection` outputs a `Mask` (of type `T <: AbstractContainer`)
containing a random instance of said type in the given `container`.

    RandomSelection{T}() where {T <: AbstractContainer}
    
The state mode of RandomSelection `M` is forced to be Stateless.

# Examples
```jldoctest
julia> sele = RandomSelection{Residue}()
RandomSelection{ProtoSyn.Stateless,Residue}()
```
"""
struct RandomSelection{M, T} <: AbstractSelection

    sele::Opt{AbstractSelection}
    
    RandomSelection{T}(sele::Opt{AbstractSelection}) where {T <: AbstractContainer} = begin
        new{Stateless, T}(sele)
    end
end

RandomSelection{T}() where {T <: AbstractContainer} = begin
    RandomSelection{T}(nothing)
end

state_mode_type(::RandomSelection{M, T}) where {M, T} = M
selection_type(::RandomSelection{M, T}) where {M, T} = T

# --- RandomSelection Select ---------------------------------------------------
function select(sele::RandomSelection{Stateless, T}, container::AbstractContainer) where {T <: AbstractContainer}

    if sele.sele === nothing
        n_items               = counter(T)(container)
        mask                  = Mask{T}(n_items)
        mask[rand(1:n_items)] = true
        return mask
    else
        n_items = counter(T)(container)
        mask = Mask{T}(n_items)
        _selector = sele.sele(container)
        mask[rand(findall(_selector.content))] = true
        return mask
    end
end

# --- Short Syntax -------------------------------------------------------------
export @ran_str
macro ran_str(p, flags...); RandomSelection{Atom}(FieldSelection{Atom}(p, :name; is_regex = parse_flags(flags))); end


# ------------------------------------------------------------------------------
export RandomRangeSelection
# Note: RandomRangeSelection is a LEAF selection.

"""
    RandomRangeSelection{M, T} <: AbstractSelection

A `RandomSelection` outputs a `Mask` (of type `T <: AbstractContainer`)
containing a random range of instances of said type in the given `container`.
The considered range is based on the instance's `:id`.

    RandomRangeSelection{T}() where {T <: AbstractContainer}
    
The state mode of RandomRangeSelection `M` is forced to be Stateless.

# Examples
```jldoctest
julia> sele = RandomRangeSelection{Residue}()
RandomRangeSelection{ProtoSyn.Stateless,Residue}()
```
"""
struct RandomRangeSelection{M, T} <: AbstractSelection
    
    RandomRangeSelection{T}() where {T <: AbstractContainer} = begin
        new{Stateless, T}()
    end
end

state_mode_type(::RandomRangeSelection{M, T}) where {M, T} = M
selection_type(::RandomRangeSelection{M, T}) where {M, T} = T

# --- RandomSelection Select ---------------------------------------------------
function select(::RandomRangeSelection{Stateless, T}, container::AbstractContainer) where {T <: AbstractContainer}

    n_items         = counter(T)(container)
    mask            = Mask{T}(n_items)
    a               = rand(1:n_items)
    b               = rand(1:n_items)
    _min, _max      = min(a, b), max(a, b)
    mask[_min:_max] = true
    return mask
end