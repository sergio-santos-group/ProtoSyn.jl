export TrueSelection
mutable struct TrueSelection{M, T} <: AbstractSelection
    is_exit_node::Bool

    TrueSelection{T}() where {T <: AbstractContainer} = begin
        new{Stateless, T}(true)
    end
end

select(::TrueSelection{Stateless, T}, container::AbstractContainer) where {T <: AbstractContainer} = Mask{T}(trues(counter(T)(container)))

# --- Examples -----------------------------------------------------------------
# sele = TrueSelection{Atom}()