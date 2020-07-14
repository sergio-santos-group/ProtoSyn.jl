export FieldSelection
mutable struct FieldSelection{S, T} <: AbstractSelection
    is_exit_node::Bool
    pattern::Regex
    field::Symbol
    
    FieldSelection{T}(pattern::String, field::Symbol) where {T <: AbstractContainer} = begin
        new{Stateless, T}(true, Regex(pattern), field)
    end
end

state_mode_type(::FieldSelection{M, T}) where {M, T} = M
selection_type(::FieldSelection{M, T}) where {M, T} = T

# --- Select -------------------------------------------------------------------
function select(sele::FieldSelection{Stateless, T}, container::AbstractContainer) where {T <: AbstractContainer}

    n_items = counter(T)(container)
    mask = Mask{T}(n_items)

    for item in iterator(T)(container)
        if occursin(sele.pattern, getproperty(item, sele.field))
            mask[item.index] = true
        end
    end
    return mask
end

# --- Short Syntax -------------------------------------------------------------

export @sn_str, @rn_str, @an_str, @as_str
macro sn_str(p); FieldSelection{Segment}(p, :name); end
macro rn_str(p); FieldSelection{Residue}(p, :name); end
macro an_str(p); FieldSelection{Atom}(p, :name); end
macro as_str(p); FieldSelection{Atom}(p, :symbol); end