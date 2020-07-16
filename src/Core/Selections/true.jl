export TrueSelection
# Note: TrueSelection is a LEAF selection.

"""
    TrueSelection{M, T} <: AbstractSelection

A `TrueSelection` returns a `Mask{T}` with all entries set to true.

# Examples
```jldoctest
julia> sele = TrueSelection{Atom}()
TrueSelection{ProtoSyn.Stateless,Atom}(true)
```
"""
mutable struct TrueSelection{M, T} <: AbstractSelection
    is_exit_node::Bool

    TrueSelection{T}() where {T <: AbstractContainer} = begin
        new{Stateless, T}(true)
    end
end

state_mode_type(::TrueSelection{M, T}) where {M, T} = M
selection_type(::TrueSelection{M, T})  where {M, T} = T


# --- Select -------------------------------------------------------------------
select(::TrueSelection{Stateless, T}, container::AbstractContainer) where {T <: AbstractContainer} = Mask{T}(trues(counter(T)(container)))
