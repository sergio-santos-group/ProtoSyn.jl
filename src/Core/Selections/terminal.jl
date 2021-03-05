export TerminalSelection
# Note: TerminalSelection is a LEAF selection.

"""
    TerminalSelection{M, T} <: AbstractSelection

A `TerminalSelection` returns a `Mask{Residue}` selecting only the terminal
residues in a Pose/Abstract Container. Terminal residues are considered when
either:

- Are children of the respective pose's origin residue;
- Are residues without children;

# Examples
```jldoctest
julia> sele = TerminalSelection{Residue}()
TerminalSelection{ProtoSyn.Stateless,Residue}()
```
"""
mutable struct TerminalSelection{M, Residue} <: AbstractSelection
    TerminalSelection{Residue}() = begin
        new{Stateless, Residue}()
    end
end

state_mode_type(::TerminalSelection{M, T}) where {M, T} = M
selection_type(::TerminalSelection{M, T})  where {M, T} = T


# --- Select -------------------------------------------------------------------
# select(::TrueSelection{Stateless, T}, container::AbstractContainer) where {T <: AbstractContainer} = Mask{T}(trues(counter(T)(container)))

function select(sele::TerminalSelection{Stateless, Residue}, container::AbstractContainer)
    @assert typeof(container) > Residue "Can't apply TerminalSelection{Residue} to container of type $(typeof(container))"
    
    origin = ProtoSyn.origin(container).container
    n_residues = counter(Residue)(container)
    mask = Mask{Residue}(n_residues)

    for res in iterator(Residue)(container)
        if res.parent == origin || length(res.children) == 0
            mask[res.index] = true
        end
    end
    return mask
end