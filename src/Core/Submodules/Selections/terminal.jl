export TerminalSelection
# Note: TerminalSelection is a LEAF selection.

"""
    TerminalSelection()

A [`TerminalSelection`](@ref) returns a [`Mask`](@ref) selecting only the
terminal [`Residue`](@ref) instances in a [`Pose`](@ref) or `AbstractContainer`.
Terminal [`Residue`](@ref) instancs are considered when either:

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
    
    _root = ProtoSyn.root(container).container
    n_residues = counter(Residue)(container)
    mask = Mask{Residue}(n_residues)

    for res in iterator(Residue)(container)
        if res.parent == _root || length(res.children) == 0
            mask[res.index] = true
        end
    end
    return mask
end

# --- Show ---------------------------------------------------------------------

Base.show(io::IO, ts::TerminalSelection) = begin
    ProtoSyn.show(io, ts)
end

function show(io::IO, ts::TerminalSelection{M, T}, level_code::Opt{LevelCode} = nothing) where {M, T}
    lead = ProtoSyn.get_lead(level_code)
    if level_code === nothing
        level_code = LevelCode()
    end
    println(io, lead*"TerminalSelection ($T)")
end