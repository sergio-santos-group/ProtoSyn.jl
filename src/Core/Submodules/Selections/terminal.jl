export UpstreamTerminalSelection
# Note: UpstreamTerminalSelection is a LEAF selection.

"""
    UpstreamTerminalSelection{T}() where {T <: AbstractContainer}

An [`UpstreamTerminalSelection`](@ref) returns a [`Mask`](@ref) selecting only
the upstream terminal [`Residue`](@ref) or [`Atom`](@ref) instances in a
[`Pose`](@ref) or `AbstractContainer`. Upstream terminal instances are defined
as being children of the [`Pose`](@ref) or `AbstractContainer`
[`root`](@ref ProtoSyn.root).

# State mode

The state mode of [`UpstreamTerminalSelection`](@ref) `M` is forced to be `Stateless`.

# Selection type

The selection type of [`UpstreamTerminalSelection`](@ref) `T` is forced to be
either [`Residue`](@ref) or [`Atom`](@ref).

!!! ukw "Note:"
    This selection does not have a short syntax version.

# Examples
```jldoctest
julia> sele = UpstreamTerminalSelection{Residue}()
UpstreamTerminalSelection (Residue)
```
"""
mutable struct UpstreamTerminalSelection{M, T} <: AbstractSelection
    UpstreamTerminalSelection{T}() where {T <: AbstractContainer} = begin
        @assert T in [Atom, Residue] "UpstreamTerminalSelection can only select either Residue or Atom instances."
        new{Stateless, T}()
    end
end

state_mode_type(::UpstreamTerminalSelection{M, T}) where {M, T} = M
selection_type(::UpstreamTerminalSelection{M, T})  where {M, T} = T


# --- Select -------------------------------------------------------------------

function select(sele::UpstreamTerminalSelection{Stateless, T}, container::AbstractContainer) where {T <: AbstractContainer}
    @assert typeof(container) >= Residue "Can't apply UpstreamTerminalSelection{Residue} to container of type $(typeof(container))"
    
    _root       = ProtoSyn.root(container)
    _root       = T === Residue ? ProtoSyn.root(container).container : _root
    n_instances = counter(T)(container)
    mask        = Mask{T}(n_instances)

    iter = typeof(container) > T ? iterator(T)(container) : [container]
    for (index, int) in enumerate(iter)
        if int.parent == _root
            mask[index] = true
        end
    end

    return mask
end

# --- Show ---------------------------------------------------------------------
function Base.show(io::IO, uts::UpstreamTerminalSelection{M, T}, level_code::Opt{LevelCode} = nothing) where {M, T}
    lead = ProtoSyn.get_lead(level_code)
    if level_code === nothing
        level_code = LevelCode()
    end
    println(io, lead*"UpstreamTerminalSelection ($T)")
end

# ------------------------------------------------------------------------------

export DownstreamTerminalSelection
# Note: DownstreamTerminalSelection is a LEAF selection.

"""
    DownstreamTerminalSelection{T}() where {T <: AbstractContainer}

A [`DownstreamTerminalSelection`](@ref) returns a [`Mask`](@ref) selecting only
downstream terminal [`Residue`](@ref) or [`Atom`](@ref) instances in a
[`Pose`](@ref) or `AbstractContainer`. Downstream terminal instances are defined
as instances with a parent but no children.

# State mode

The state mode of [`DownstreamTerminalSelection`](@ref) `M` is forced to be `Stateless`.

# Selection type

The selection type of [`DownstreamTerminalSelection`](@ref) `T` is forced to be either
[`Residue`](@ref) or [`Atom`](@ref).

!!! ukw "Note:"
    This selection does not have a short syntax version.

# Examples
```jldoctest
julia> sele = DownstreamTerminalSelection()
DownstreamTerminalSelection (Residue)
```
"""
mutable struct DownstreamTerminalSelection{M, T} <: AbstractSelection
    DownstreamTerminalSelection{T}() where {T <: AbstractContainer} = begin
        @assert T in [Atom, Residue] "DownstreamTerminalSelection can only select either Residue or Atom instances."
        new{Stateless, T}()
    end
end

state_mode_type(::DownstreamTerminalSelection{M, T}) where {M, T} = M
selection_type(::DownstreamTerminalSelection{M, T})  where {M, T} = T


# --- Select -------------------------------------------------------------------

function select(sele::DownstreamTerminalSelection{Stateless, T}, container::AbstractContainer) where {T <: AbstractContainer}
    @assert typeof(container) >= Residue "Can't apply DownstreamTerminalSelection{Residue}() to container of type $(typeof(container))"
    
    n_instances = counter(T)(container)
    mask        = Mask{T}(n_instances)

    iter = typeof(container) > T ? iterator(T)(container) : [container]
    for (index, int) in enumerate(iter)
        if hasparent(int) && length(int.children) == 0
            mask[index] = true
        end
    end

    return mask
end

# --- Show ---------------------------------------------------------------------
function Base.show(io::IO, ts::DownstreamTerminalSelection{M, T}, level_code::Opt{LevelCode} = nothing) where {M, T}
    lead = ProtoSyn.get_lead(level_code)
    if level_code === nothing
        level_code = LevelCode()
    end
    println(io, lead*"DownstreamTerminalSelection ($T)")
end