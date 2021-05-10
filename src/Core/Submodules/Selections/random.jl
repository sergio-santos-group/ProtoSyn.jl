export RandomSelection
# Note: RandomSelection is a LEAF selection.

"""
    RandomSelection{T}() where {T <: AbstractContainer}

A [`RandomSelection`](@ref) outputs a [`Mask`](@ref) (of type
`T <: AbstractContainer`) containing a random instance of said type in the given
`container`.

# State mode    

The state mode of [`RandomSelection`](@ref) `M` is forced to be `Stateless`.

# Selection type

The selection type of [`RandomSelection`](@ref) can be any
`T <: AbstractContainer`.

!!! ukw "Note:"
    This selection does not have a short syntax version.

# Examples
```jldoctest
julia> sele = RandomSelection{Residue}()
RandomSelection › Residue.id
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

# --- Show ---------------------------------------------------------------------

Base.show(io::IO, rs::RandomSelection) = begin
    ProtoSyn.show(io, rs)
end

function show(io::IO, rs::RandomSelection{M, T}, levels::Opt{BitArray} = nothing) where {M, T}
    lead = ProtoSyn.get_lead(levels)
    if levels === nothing
        levels = BitArray([])
    end
    println(io, lead*"RandomSelection › $T.id")
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
    
The state mode of [`RandomRangeSelection`](@ref) `M` is forced to be `Stateless`.

# Selection type

The selection type of [`RandomRangeSelection`](@ref) can be any
`T <: AbstractContainer`.

!!! ukw "Note:"
    This selection does not have a short syntax version.

# Examples
```jldoctest
julia> sele = RandomRangeSelection{Residue}()
RandomRangeSelection › Residue.id
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

# --- Show ---------------------------------------------------------------------

Base.show(io::IO, rrs::RandomRangeSelection) = begin
    ProtoSyn.show(io, rrs)
end

function show(io::IO, rrs::RandomRangeSelection{M, T}, levels::Opt{BitArray} = nothing) where {M, T}
    lead = ProtoSyn.get_lead(levels)
    if levels === nothing
        levels = BitArray([])
    end
    println(io, lead*"RandomRangeSelection › $T.id")
end

# ------------------------------------------------------------------------------
export RandomSelectionFromList
# Note: RandomSelectionFromList is a BRANCH selection.

"""
    RandomSelectionFromList(selections::Vector{T}) where {T <: AbstractSelection}

A [`RandomSelectionFromList`](@ref) outputs a [`Mask`](@ref) (of type
`T <: AbstractContainer`). This [`Mask`](@ref) is the result of the application
of a randomly selected `AbstractSelection` from the provided list of
`AbstractSelection` instances `selections`.

!!! ukw "Note:"
    All the given `AbstractSelection` instances must be of the same type.

# State mode

The state mode of [`RandomSelectionFromList`](@ref) `M` is the same as the state
mode of the provided list of `AbstractSelection` instances (which are all of the
same type).

# Selection type

The selection type of [`RandomSelectionFromList`](@ref) `M` is the same as the
selection type of the provided list of `AbstractSelection` instances (which are
all of the same type).

!!! ukw "Note:"
    This selection does not have a short syntax version.

# Examples
```jldoctest
julia> s = ProtoSyn.RandomSelectionFromList([rid"1", rid"2"])
RandomSelectionFromList ❯ (Residue)
 ├── SerialSelection › Residue.id = 1
 └── SerialSelection › Residue.id = 2
 
julia> s = ProtoSyn.RandomSelectionFromList([rid"1", rn"CBZ"])
ERROR: AssertionError: RandomSelectionFromList `selections` elements must be all of the same type.
```
"""
struct RandomSelectionFromList{M, T} <: AbstractSelection
    
    selections::Vector{AbstractSelection}

    RandomSelectionFromList(selections::Vector{T}) where {T <: AbstractSelection} = begin
        @assert length(Set(map(i -> typeof(i), selections))) == 1 "RandomSelectionFromList `selections` elements must be all of the same type."
        new{state_mode_type(selections[1]), selection_type(selections[1])}(selections)
    end
end

state_mode_type(s::RandomSelectionFromList{M, T}) where {M, T} = M
selection_type(s::RandomSelectionFromList{M, T}) where {M, T} = T

# --- RandomSelectionFromList Select -------------------------------------------
function select(rsfl::RandomSelectionFromList{M, T}, container::AbstractContainer) where {M, T <: AbstractContainer}
    return select(rsfl.selections[rand(1:end)], container)
end

# --- Show ---------------------------------------------------------------------

Base.show(io::IO, rsfl::RandomSelectionFromList) = begin
    ProtoSyn.show(io, rsfl)
end

function show(io::IO, rsfl::RandomSelectionFromList{M, T}, levels::Opt{BitArray} = nothing) where {M, T}
    lead = ProtoSyn.get_lead(levels)
    if levels === nothing
        levels = BitArray([])
    end
    println(io, lead*"RandomSelectionFromList ❯ ($T)")
    for selection in rsfl.selections[1:(end-1)]
        ProtoSyn.show(io, selection, vcat(levels, true))
    end
    ProtoSyn.show(io, rsfl.selections[end], vcat(levels, false))
end