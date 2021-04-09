export RandomSelection
# Note: RandomSelection is a LEAF selection.

"""
    RandomSelection{T}() where {T <: AbstractContainer}

A [`RandomSelection`](@ref) outputs a [`Mask`](@ref) (of type
`T <: AbstractContainer`) containing a random instance of said type in the given
`container`.

# State mode    

The state mode of [`RandomSelection`] `M` is forced to be `Stateless`.

# Selection type

The selection type of [`RandomSelection`](@ref) can be any
`T <: AbstractContainer`.

!!! ukw "Note:"
    This selection does not have a short syntax version.

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
        # _selector = sele.sele(container)
        _selector = promote(sele.sele, T)(container)
        mask[rand(findall(_selector.content))] = true
        return mask
    end
end

# ------------------------------------------------------------------------------
export RandomRangeSelection
# Note: RandomRangeSelection is a LEAF selection.

"""
    RandomRangeSelection{T}() where {T <: AbstractContainer}

A [`RandomRangeSelection`](@ref) outputs a [`Mask`](@ref) (of type
`T <: AbstractContainer`) containing a random range of instances of said type in
the given `container`. The considered range is based on the instance's `:id`.

# State mode
    
The state mode of [`RandomRangeSelection`] `M` is forced to be `Stateless`.

# Selection type

The selection type of [`RandomRangeSelection`](@ref) can be any
`T <: AbstractContainer`.

!!! ukw "Note:"
    This selection does not have a short syntax version.

# Examples
```jldoctest
julia> sele = RandomRangeSelection{Residue}()
RandomRangeSelection{ProtoSyn.Stateless,Residue}()
```

!!! ukw "Note:"
    This selection assumes that all `Abstractcontainer` instances are ordered
    (i.e: a random range between atom 1 and atom 10 will select atoms 1, 2, 3,
    4, 5, 6, 7, 8, 9 and 10).

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