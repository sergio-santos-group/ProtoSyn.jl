export TerminalSelection
# Note: TerminalSelection is a LEAF selection.

"""
    TerminalSelection()

A [`TerminalSelection`](@ref) returns a [`Mask`](@ref) selecting only the
terminal [`Residue`](@ref) instances in a [`Pose`](@ref) or `AbstractContainer`.
Terminal [`Residue`](@ref) instances are considered when either:

- Are children of the respective pose's root residue;
- Are residues without children;

# State mode

The state mode of [`TerminalSelection`](@ref) `M` is forced to be `Stateless`.

# Selection type

The selection type of [`TerminalSelection`](@ref) `T` is forced to be [`Residue`](@ref).

!!! ukw "Note:"
    This selection does not have a short syntax version.

# Examples
```jldoctest
julia> sele = TerminalSelection()
TerminalSelection (Residue)
```
"""
mutable struct TerminalSelection{M, T} <: AbstractSelection
    TerminalSelection() = begin
        new{Stateless, Residue}()
    end
end

state_mode_type(::TerminalSelection{M, T}) where {M, T} = M
selection_type(::TerminalSelection{M, T})  where {M, T} = T


# --- Select -------------------------------------------------------------------

function select(sele::TerminalSelection{Stateless, Residue}, container::AbstractContainer)
    @assert typeof(container) > Residue "Can't apply TerminalSelection{Residue} to container of type $(typeof(container))"
    
    _root      = ProtoSyn.root(container).container
    n_residues = counter(Residue)(container)
    mask       = Mask{Residue}(n_residues)

    for res in iterator(Residue)(container)
        if res.parent == _root || length(res.children) == 0
            mask[res.index] = true
        end
    end
    return mask
end

# --- Show ---------------------------------------------------------------------
function Base.show(io::IO, ts::TerminalSelection{M, T}, level_code::Opt{LevelCode} = nothing) where {M, T}
    lead = ProtoSyn.get_lead(level_code)
    if level_code === nothing
        level_code = LevelCode()
    end
    println(io, lead*"TerminalSelection ($T)")
end

# ------------------------------------------------------------------------------

export NTerminalSelection
# Note: NTerminalSelection is a LEAF selection.

"""
    NTerminalSelection()

A [`NTerminalSelection`](@ref) returns a [`Mask`](@ref) selecting only the N
terminal [`Residue`](@ref) instances in a [`Pose`](@ref) or `AbstractContainer`.
N terminal [`Residue`](@ref) instances are considered when they are residues
without children;

# State mode

The state mode of [`NTerminalSelection`](@ref) `M` is forced to be `Stateless`.

# Selection type

The selection type of [`NTerminalSelection`](@ref) `T` is forced to be
[`Residue`](@ref).

!!! ukw "Note:"
    This selection does not have a short syntax version.

# Examples
```jldoctest
julia> sele = NTerminalSelection()
N-TerminalSelection (Residue)
```
"""
mutable struct NTerminalSelection{M, T} <: AbstractSelection
    NTerminalSelection() = begin
        new{Stateless, Residue}()
    end
end

state_mode_type(::NTerminalSelection{M, T}) where {M, T} = M
selection_type(::NTerminalSelection{M, T})  where {M, T} = T


# --- Select -------------------------------------------------------------------

function select(sele::NTerminalSelection{Stateless, Residue}, container::AbstractContainer)
    @assert typeof(container) > Residue "Can't apply NTerminalSelection{Residue} to container of type $(typeof(container))"
    
    _root      = ProtoSyn.root(container).container
    n_residues = counter(Residue)(container)
    mask       = Mask{Residue}(n_residues)

    for (index, res) in enumerate(iterator(Residue)(container))
        if res.parent == _root
            mask[index] = true
        end
    end
    return mask
end

# --- Show ---------------------------------------------------------------------
function Base.show(io::IO, ts::NTerminalSelection{M, T}, level_code::Opt{LevelCode} = nothing) where {M, T}
    lead = ProtoSyn.get_lead(level_code)
    if level_code === nothing
        level_code = LevelCode()
    end
    println(io, lead*"N-TerminalSelection ($T)")
end

# ------------------------------------------------------------------------------

export CTerminalSelection
# Note: CTerminalSelection is a LEAF selection.

"""
    CTerminalSelection()

A [`CTerminalSelection`](@ref) returns a [`Mask`](@ref) selecting only the C
terminal [`Residue`](@ref) instances in a [`Pose`](@ref) or `AbstractContainer`.
C terminal [`Residue`](@ref) instances are considered when they are children of
the respective pose's [`root`](@ref ProtoSyn.root) residue;

# State mode

The state mode of [`CTerminalSelection`](@ref) `M` is forced to be `Stateless`.

# Selection type

The selection type of [`CTerminalSelection`](@ref) `T` is forced to be
[`Residue`](@ref).

!!! ukw "Note:"
    This selection does not have a short syntax version.

# Examples
```jldoctest
julia> sele = CTerminalSelection()
C-TerminalSelection (Residue)
```
"""
mutable struct CTerminalSelection{M, T} <: AbstractSelection
    CTerminalSelection() = begin
        new{Stateless, Residue}()
    end
end

state_mode_type(::CTerminalSelection{M, T}) where {M, T} = M
selection_type(::CTerminalSelection{M, T})  where {M, T} = T


# --- Select -------------------------------------------------------------------

function select(sele::CTerminalSelection{Stateless, Residue}, container::AbstractContainer)
    @assert typeof(container) > Residue "Can't apply CTerminalSelection{Residue} to container of type $(typeof(container))"
    
    _root      = ProtoSyn.root(container).container
    n_residues = counter(Residue)(container)
    mask       = Mask{Residue}(n_residues)

    for (index, res) in enumerate(iterator(Residue)(container))
        if length(res.children) == 0
            mask[index] = true
        end
    end
    return mask
end

# --- Show ---------------------------------------------------------------------
function Base.show(io::IO, ts::CTerminalSelection{M, T}, level_code::Opt{LevelCode} = nothing) where {M, T}
    lead = ProtoSyn.get_lead(level_code)
    if level_code === nothing
        level_code = LevelCode()
    end
    println(io, lead*"C-TerminalSelection ($T)")
end